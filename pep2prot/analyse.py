import numpy as np
import pandas as pd
from pathlib import Path

from .read import read_isoquant_peptide_report, read_fastas
from .graphs.pep_prot_graph import ProtPepGraph
from .preprocessing import preprocess_isoquant_peptide_report, get_protein_coverages, complex_cluster_buster, simple_cluster_buster 
from .intensities import get_prot_intensities
from .postprocessing import summarize_prots, get_stats, prettify_protein_informations, get_full_report


def isoquant_peptide_report(pep_rep_path, 
                            fastas_path,
                            min_pepNo_per_prot=2,
                            verbose=False,
                            cluster_buster=complex_cluster_buster,
                            full_outcome=False):
    """Analyze IsoQuant peptide report.

    Args:
        pep_rep_path (str or pathlib.Path): Path to a peptide report.
        fastas_path (str or pathlib.Path): Path to a fasta file.
        min_pepNo_per_prot (int): The minimal number of peptides per protein.
        verbose (boolean): Show messages?
        cluster_buster (function): The clustering procedure to use to get rid of peptides spread over different clusters.
        full_outcome (boolean): Return more output?
    Returns:
        tuple: first two elements will always be the prettified intensity reports.
        If 'full_outcome' is true, the result will contain after these:
        the full protein-peptide graph, its induced minimal graph, 
        tuple with lonely peptides and proteins,
        tuple with peptides attached only to proteins with low support and these proteins,
        proteins groups not needed to cover existing peptide groups,
        minimal, deconvoluted, and maximal intensities of protein groups.
    """
    pep_rep_path = Path(pep_rep_path).expanduser()
    fastas_path  = Path(fastas_path).expanduser()

    if verbose:
        print('Reading isoquant report.')
    D = read_isoquant_peptide_report(pep_rep_path)

    if verbose:
        print('Preprocessing isoquant report.')
    D, I_cols = preprocess_isoquant_peptide_report(D)

    if verbose:
        print('Reading fastas.')
    fastas = read_fastas(fastas_path)
    observed_prots = {r for rg in D.prots for r in rg}
    assert all(r in fastas.index for r in observed_prots), "It seems that you are using a different set of fastas than the peptide annotation software before. Repent please."

    if verbose:
        print('Calculating coverages.')
    prots = get_protein_coverages(D, fastas)
    uni_cols  = ['peptide_overall_max_score','peptide_fdr_level',
                 'peptide_overall_replication_rate','prots',
                 'pre_homology_accessions','pi','mw']
    if verbose:
        print('Aggregating the same peptides that ended up in different clusters.')
    DD = cluster_buster(D, I_cols, uni_cols) # agg same peptides in various clusters
    assert np.all(DD.groupby(['pep','prots']).size() == 1), "Some peptides are still across different clusters."

    if verbose:
        print('Building the peptide-protein graph.')

    G = ProtPepGraph((r,p) for p, rg in zip(DD.index, DD.prots) for r in rg)
    lonely, unsupported = G.remove_lonely_and_unsupported(min_pepNo_per_prot)
    H, rejected = G.get_minimal_graph()

    pep2pepgr = {p:pg for pg in H.peps() for p in pg}
    DDinH = DD.loc[pep2pepgr] # peps in H: no simple prot-pep pairs, no unnecessary prots?
    DDinH['pepgr'] = DDinH.index.map(pep2pepgr)
    pep_I = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities

    if verbose:
        print('Spreading intensities from peptides to proteins.')
    prot_min_I, prot_I, prot_max_I = get_prot_intensities(H, pep_I)
    prot_info = summarize_prots(H, fastas, prots.pep_coverage)
    prot_I_nice = prettify_protein_informations(prot_I, prot_info)
    prot_stats = get_stats(prot_min_I, prot_I, prot_max_I)
    
    if verbose:
        print('Preparing reports.')
    all_prots = get_full_report(prot_min_I, prot_I, prot_max_I)
    all_prot_nice = prettify_protein_informations(all_prots, prot_info)

    if full_outcome:
        return prot_I_nice, all_prot_nice, G, H, lonely, unsupported, rejected, prot_min_I, prot_I, prot_max_I
    else:
        return prot_I_nice, all_prot_nice

