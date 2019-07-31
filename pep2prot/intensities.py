import numpy as np
import pandas as pd


debug = True

def get_prot_intensities(H, peps_I):
    """Get intensities based on a protein-peptide graph.

    Args:
        H (ProtPepGraph): Ready protein-peptide graph.
        peps_I (pd.DataFrame): Observed intensities of individual peptides.
    Returns:
        Minimal, deconvoluted, and maximal protein intensities.
    """
    prot_pep = ('prot','pep')
    Hdf = pd.DataFrame.from_records(H.prot_pep_pairs(), columns=prot_pep, index=prot_pep)
    peps2prots_max_I = Hdf.join(peps_I, on='pep')
    peps2prots = peps2prots_max_I.index
    prots_max_I = peps2prots_max_I.groupby('prot').sum()
    if debug:
        print('Some NaNs in "prots_max_I":', peps2prots_max_I.isna().any().any())

    # peptides that point uniquely to one protein group
    unipeps = pd.Index((pg for pg in H.peps() if H.degree(pg) == 1), name='pep')
    
    # intensities going from unique peptides towards their respective proteins
    unipeps2prots_I = peps2prots_max_I[peps2prots.get_level_values('pep').isin(unipeps)]
    if debug:
        print('Some NaNs in "unipeps2prots_I":', unipeps2prots_I.isna().any().any())

    # there is one unique peptide per protein, so _I == _max_I == _min_I
    uniprots_min_I = unipeps2prots_I.reset_index('pep', drop=True)
    
    # prots with some uniquely identified peptides
    uniprots = uniprots_min_I.index
    
    # blade = razor??? probably not, so stick with blade.
    bladeprots = pd.Index((r for r in H.prots() if not r in uniprots), name='prot')
    bladeprots_zero_I = pd.DataFrame(np.zeros(shape=(len(bladeprots),
                                                     len(unipeps2prots_I.columns))),
                                     index=bladeprots, columns=unipeps2prots_I.columns)
    # prots = uniprots ⊔ bladeprots
    prots_min_I = pd.concat([uniprots_min_I, bladeprots_zero_I])
    if debug:
        print('Some NaNs in "prots_min_I":', unipeps2prots_I.isna().any().any())

    sorted_prots = np.sort(prots_min_I.index) # sorting to compare intensities
    prots_min_I = prots_min_I.loc[sorted_prots]
    prots_max_I = prots_max_I.loc[sorted_prots]
    assert np.all(prots_min_I <= prots_max_I), "Some minimal intensities are not smaller then the maximal ones. Report to Matteo."

    uniprots_curr_I = uniprots_min_I # reusable: intensity from unique peps pushed to uniprots
    bladepeps = pd.Index({p for r in bladeprots 
                            for p in H[r] if all(rr in bladeprots for rr in H[p])},
                         name='pep') # peps with prots without unique peps
    bladepeps_I = peps_I.loc[bladepeps]
    bladepeps_protcnts = pd.Series((H.degree(p) for p in bladepeps), index=bladepeps)

    # distribute blade peps intensity uniformly among neighbor prots: hence the division below
    bladepeps2bladeprots_I = pd.DataFrame.from_records(
        ((r,p) for p in bladepeps for r in H[p]), 
        columns=prot_pep,
        index=prot_pep).join( bladepeps_I.div(bladepeps_protcnts, axis='index'), on='pep')
    bladeprots_curr_I = bladepeps2bladeprots_I.groupby('prot').sum()

    # bladeprots might receive some more intensity from other peptides, hence _curr_I.
    # peps = unipeps ⊔ bladepeps ⊔ otherpeps
    # unipeps neighbor uniprots, bladepeps neighbor bladeprots, otherpeps neighbor some prots from both sets
    otherpeps = pd.Index({p for p in H.peps() if p not in unipeps and p not in bladepeps}, name='pep')
    otherpeps_I = peps_I.loc[otherpeps]

    # above intensities will be distributed proportionally to intensities prots received from bladepeps and unipeps.
    # call these the current intensities, curr_I.
    # bladeprots and uniprots are disjoint sets, so we concat them
    prots_curr_I = pd.concat([uniprots_curr_I, bladeprots_curr_I])
    
    # populating edges with other intensities from otherpeps. Let's call their prots mixprots
    otherpeps2mixprots_I = pd.DataFrame(index=peps2prots).join(otherpeps_I, on='pep', how='right')
    otherpeps2mixprots = otherpeps2mixprots_I.index
    mixprots = otherpeps2mixprots_I.index.get_level_values('prot').unique()

    # need weights for otherpeps2mixprots_I
    weights = pd.DataFrame(index=otherpeps2mixprots).join(prots_curr_I, on='prot')
    eps = np.finfo(type(weights.iloc[0,0])).eps# machine precision
    weights += eps # whenever we have 0,0,0, we spread intensities proportionally to eps/3eps = 1/3, rather than 0/0.
    weights_mixprot_I = weights.groupby('pep').sum()
    weights = weights.div(weights_mixprot_I, axis='index')

    assert not np.any(weights.isnull()), "Weight cannot result in any NaNs."
    otherpeps2mixprots_I = otherpeps2mixprots_I * weights
    prots_curr_I.loc[mixprots] += otherpeps2mixprots_I.groupby('prot').sum()# update only mixprots
    prots_I = prots_curr_I.loc[sorted_prots]

    assert np.all(prots_min_I <= prots_I), "Some deconvoluted intensities are smaller then minimal intensities."
    assert np.all(prots_I <= prots_max_I), "Some deconvoluted intensities are larger then maximal intensities."

    return prots_min_I, prots_I, prots_max_I

