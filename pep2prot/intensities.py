import numpy as np
import pandas as pd
from pandas import DataFrame as df


debug = False

def get_prot_intensities(H, pep_I):
    """Get intensities based on a protein-peptide graph.

    Args:
        H (ProtPepGraph): Ready protein-peptide graph.
        pep_I (pd.DataFrame): Observed intensities of individual peptides.
    Returns:
        Minimal, deconvoluted, and maximal protein intensities.
    """
    (prot_cnt, pep_cnt), pep2prot_cnt = H.nodes_cnt(), len(H.edges)
    _, Icols_cnt = pep_I.shape
    Icols = pep_I.columns
    assert pep_cnt == _, "H has {} peptides and pep_I {} peptides!".format(pep_cnt, _)

    pep2prot = pd.MultiIndex.from_tuples(H.prot_pep_pairs(), names=('prot','pep'))
    pep2prot_maxI = df(index=pep2prot).join(pep_I, on='pep')
    prot_maxI = pep2prot_maxI.groupby('prot').sum()
    
    assert pep2prot_maxI.shape[0] == pep2prot_cnt, "Row number in pep2prot_maxI is not equal to the number of edges in the graph."
    assert prot_maxI.shape[0] == prot_cnt, "Row number in prot_maxI is not equal to the number of protein groups."
    assert not np.any(prot_maxI.isna()), "Missing values in the maximal intensities."
    
    # peptides unique for one protein group
    unipep = pd.Index(H.peps(deg=1), name='pep')

    # intensities going from unique peptides towards their respective proteins
    unipep2prot_I = pep2prot_maxI[pep2prot.get_level_values('pep').isin(unipep)]    
    uniprot_minI = unipep2prot_I.groupby('prot').sum()

    # proteins with at least one uniquely identifying peptide group
    uniprot = uniprot_minI.index
    
    # proteins that have no uniquely identifying peptides (uni backwards is inu)
    inuprot = pd.Index((r for r in H.prots() if not r in uniprot), name='prot')
    inuprot_noI = df(np.zeros(shape=(len(inuprot), Icols_cnt)), index=inuprot, columns=Icols)

    # prot = uniprot ⊔ inuprot
    prot_minI = pd.concat([uniprot_minI, inuprot_noI])
    assert prot_minI.shape == prot_maxI.shape, "Incompatible shapes of minimal and maximal intensities."

    sorted_prot = np.sort(prot_minI.index) # sorting to compare intensities
    prot_minI = prot_minI.loc[sorted_prot]
    prot_maxI = prot_maxI.loc[sorted_prot]
    assert np.all(prot_minI <= prot_maxI), "Some minimal intensities are not smaller then the maximal ones. Report to Matteo."



    ####### deconvoluting intensities #######

    # peptides surrounded by inuprots, i.e. their proteins are all not directly identified by some peptide.
    inupep = pd.Index({p for r in inuprot for p in H[r] if all(R in inuprot for R in H[p])}, name='pep')
    inupep_I = pep_I.loc[inupep]
    inupep_protcnts = [H.degree(p) for p in inupep]
    assert all(cnt > 0 for cnt in inupep_protcnts)
    # intensity is divided uniformly among the inuprots
    inupep_divI = inupep_I.div(inupep_protcnts, axis='index')
    inupep2inuprot_I = df.from_records(((r,p) for p in inupep for r in H[p]), 
                                       columns=('prot','pep'), index=('prot','pep'))
    inupep2inuprot_I = inupep2inuprot_I.join(inupep_divI, on='pep')
    inuprot_I = inupep2inuprot_I.groupby('prot').sum()

    # pep = unipep ⊔ inupep ⊔ bpep (from bridge pep: they bridge uniprots and inuprots)
    # unipep neighbors uniprot, inupep neighbors inuprot, bpep neighbors some prots from both sets
    bpep = pd.Index({p for p in H.peps() if p not in unipep and p not in inupep}, name='pep')
    bpep_I = pep_I.loc[bpep]
    prot_I = pd.concat([uniprot_minI, inuprot_I])
    assert prot_I.shape[0] == prot_cnt, "prot_I has wrong number of proteins."

    # populating edges with other intensities from bpep.
    bpep2bprot_I = df(index=pep2prot).join(bpep_I, on='pep', how='right')
    bpep2bprot = bpep2bprot_I.index

    # proteins that neighbor bpeps: some are in uniprot, other in inuprot.
    bprot = bpep2bprot.get_level_values('prot').unique()

    # need weights for bpep2bprot_I (proportionality to the already attributed intensities)
    weights = df(index=bpep2bprot).join(prot_I, on='prot')
    eps = np.finfo(type(weights.iloc[0,0])).eps# machine precision
    # avoiding division by zero for zero peptide levels in some columnss.
    # if only zeros (say, three zeros), spread I at eps/3eps = 1/3, rather than 0/0.
    weights += eps
    weights = weights.div(weights.groupby('pep').sum(), axis='index')
    assert not np.any(weights.isnull()), "Weight cannot result in any NaNs."

    bpep2bprot_I = bpep2bprot_I * weights
    bprot_I = bpep2bprot_I.groupby('prot').sum()
    assert len(bprot) == bprot_I.shape[0], "bprot_I has wrong number of proteins in rows."

    prot_I.loc[bprot] += bprot_I# update only bprot
    prot_I = prot_I.loc[sorted_prot]
    assert prot_I.shape[0] == prot_cnt, "prot_I has a wrong number of proteins."
    assert np.all(prot_minI <= prot_I), "Some deconvoluted intensities are smaller then minimal intensities."
    assert np.all(prot_I <= prot_maxI), "Some deconvoluted intensities are larger then maximal intensities."

    return prot_minI, prot_I, prot_maxI

