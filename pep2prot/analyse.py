import numpy as np
import pandas as pd
from pathlib import Path

from .read import read_isoquant_peptide_report, read_fastas
from .graphs import ProtPepGraph, BiGraph, get_peptide_protein_graph
from .preprocessing import preprocess_isoquant_peptide_report, get_protein_coverages, complex_cluster_buster, simple_cluster_buster 
from .intensities import get_prot_intensities
from .postprocessing import summarize_prots, get_stats, prettify_protein_informations, get_full_report


def isoquant_peptide_report(pep_rep_path, 
                            fastas_path,
                            cluster_buster=complex_cluster_buster,
                            verbose=False):
    """Analyze IsoQuant peptide report.

    Args:
        pep_rep_path (str or pathlib.Path): Path to a peptide report.
        fastas_path (str or pathlib.Path): Path to a fasta file.

    """
    pep_rep_path = Path(pep_rep_path).expanduser()
    fastas_path  = Path(fastas_path).expanduser()

    if verbose:
        print('Reading isoquant report.')
    D         = read_isoquant_peptide_report(pep_rep_path)

    if verbose:
        print('Preprocessing isoquant report.')
    D, I_cols = preprocess_isoquant_peptide_report(D)

    if verbose:
        print('Reading fastas.')
    fastas    = read_fastas(fastas_path, {r for rg in D.prots for r in rg})

    if verbose:
        print('Calculating coverages.')
    prots     = get_protein_coverages(D, fastas)
    uni_cols  = ['peptide_overall_max_score','peptide_fdr_level',
                      'peptide_overall_replication_rate','prots',
                      'pre_homology_accessions','pi','mw']
    if verbose:
        print('Aggregating the same peptides that ended up in different clusters.')
    DD = cluster_buster(D, I_cols, uni_cols) # agg same peptides in various clusters
    assert np.all(DD.groupby(['pep','prots']).size() == 1), "Some peptides are still across different clusters."

    if verbose:
        print('Building the peptide-protein graph.')
    H, prots_no_peps, peps_no_prots, beckham_prots = get_peptide_protein_graph(DD) 
    pep2pepgr      = {p:pg for pg in H.peps() for p in pg}
    DDinH          = DD.loc[pep2pepgr] # peps in H: no simple prot-pep pairs, no unnecessary prots?
    DDinH['pepgr'] = DDinH.index.map(pep2pepgr)
    peps_I         = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities

    if verbose:
        print('Spreading intensities from peptides to proteins.')
    prots_min_I, prots_I, prots_max_I = get_prot_intensities(H, peps_I)
    prot_info      = summarize_prots(H, fastas, prots.pep_coverage)
    prots_I_nice   = prettify_protein_informations(prots_I, prot_info)
    prots_stats    = get_stats(prots_min_I, prots_I, prots_max_I)
    
    if verbose:
        print('Preparing reports.')
    all_prots      = get_full_report(prots_min_I, prots_I, prots_max_I)
    all_prots_nice = prettify_protein_informations(all_prots, prot_info)

    return prots_I_nice, all_prots_nice

