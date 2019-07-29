%load_ext autoreload
%autoreload 2

import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter

from pep2prot.read import read_isoquant_peptide_report, read_fastas
from pep2prot.graphs import ProtPepGraph, BiGraph, get_peptide_protein_graph
from pep2prot.preprocessing import preprocess_isoquant_peptide_report, get_protein_coverages, complex_cluster_buster, simple_cluster_buster 
from pep2prot.intensities import get_prot_intensities
from pep2prot.postprocessing import summarize_prots, get_stats, prettify_protein_informations, get_full_report


test_data = Path(r"/home/matteo/Projects/pep2prot/pep2prot/data")
pep_rep_path = test_data/'peptide_report.csv'
fastas_path = test_data/'mouse.fasta'
cluster_buster=complex_cluster_buster

D = read_isoquant_peptide_report(pep_rep_path)
D, I_cols = preprocess_isoquant_peptide_report(D)
fastas = read_fastas(fastas_path)
observed_prots = {r for rg in D.prots for r in rg}
assert all(r in fastas.index for r in observed_prots), "It seems that you are using a different set of fastas than the peptide annotation software before. Repent please."
prots     = get_protein_coverages(D, fastas)
uni_cols  = ['peptide_overall_max_score','peptide_fdr_level',
                  'peptide_overall_replication_rate','prots',
                  'pre_homology_accessions','pi','mw']

DD = cluster_buster(D, I_cols, uni_cols) # agg same peptides in various clusters
assert np.all(DD.groupby(['pep','prots']).size() == 1), "Some peptides are still across different clusters."

H, prots_no_peps, peps_no_prots, beckham_prots = get_peptide_protein_graph(DD) 
Counter(len(r) for r in H.prots())
Counter(D.prots.map(len))
all_prots = {r for rg in D.prots for r in rg}
prots_in_prot_groups = {r for rg in H.prots() for r in rg}
len(prots_in_prot_groups) + len(prots_no_peps) + len(beckham_prots)
len(all_prots)

# pep2pepgr      = {p:pg for pg in H.peps() for p in pg}
# DDinH          = DD.loc[pep2pepgr] # peps in H: no simple prot-pep pairs, no unnecessary prots?
# DDinH['pepgr'] = DDinH.index.map(pep2pepgr)
# peps_I         = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities

# if verbose:
#     print('Spreading intensities from peptides to proteins.')
# prots_min_I, prots_I, prots_max_I = get_prot_intensities(H, peps_I)
# prot_info      = summarize_prots(H, fastas, prots.pep_coverage)
# prots_I_nice   = prettify_protein_informations(prots_I, prot_info)
# prots_stats    = get_stats(prots_min_I, prots_I, prots_max_I)

# if verbose:
#     print('Preparing reports.')
# all_prots      = get_full_report(prots_min_I, prots_I, prots_max_I)
# all_prots_nice = prettify_protein_informations(all_prots, prot_info)

