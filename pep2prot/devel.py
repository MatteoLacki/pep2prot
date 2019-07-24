%load_ext autoreload
%autoreload 2

from collections import Counter
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 30)#display whole column without truncation
from pathlib import Path

from pep2prot.read import read_isoquant_peptide_report, read_fastas
from pep2prot.graphs import ProtPepGraph, BiGraph, get_peptide_protein_graph
from pep2prot.preprocessing import preprocess_isoquant_peptide_report, get_protein_coverages, complex_cluster_buster, simple_cluster_buster 
from pep2prot.intensities import get_prot_intensities
from pep2prot.postprocessing import get_info_on_prots, get_stats


min_pepNo_per_prot = 2
path = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()
cluster_buster = complex_cluster_buster # simple_cluster_buster

D = read_isoquant_peptide_report(path/'peptide_report.csv')
D, I_cols = preprocess_isoquant_peptide_report(D)
unique_columns = ['peptide_overall_max_score','peptide_fdr_level',
                  'peptide_overall_replication_rate','prots',
                  'pre_homology_accessions','pi','mw']

fastas = read_fastas(path/'mouse.fasta', {r for rg in D.prots for r in rg})
prots  = get_protein_coverages(D, fastas)
# plt.hist(prots.pep_coverage, bins=100);plt.show()

DD = cluster_buster(D, I_cols, unique_columns)# agg same peptides in various clusters
assert np.all(DD.groupby(['pep','prots']).size() == 1), "Some peptides are still across different clusters."
H, prots_no_peps, peps_no_prots, beckhams_prots = get_peptide_protein_graph(DD) 

pep2pepgr = {p:pg for pg in H.peps() for p in pg}
DDinH = DD.loc[pep2pepgr] # peps in H: no simple prot-pep pairs, no unnecessary prots?
DDinH['pepgr'] = DDinH.index.map(pep2pepgr)
peps_I = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities


# pep2pepgr = {p:pg for pg in H.peps() for p in pg}
# DDinH = DD.loc[pep2pepgr]
# protgroups = set(H.prots())
# DDinH = DDinH[DDinH.prots.isin(protgroups)] # why getting rid of some groups hurts????
# peps_II = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities
#ERRROORROROROR this has to be fixed: I have no idea what's wrong here.
peps_I
# peps_II

prots_min_I, prots_I, prots_max_I = get_prot_intensities(H, peps_I)
# prots_stats = get_stats(prots_min_I, prots_I, prots_max_I)
prots_HF = prots.loc[(r for rg in prots_I.index for r in rg)]
# plt.hist(prots_HF.pep_coverage, bins=100); plt.show()

prot_info = get_info_on_prots(prots_I.index, fastas)# find a better name for that function for God's sake...
prot2pep            = pd.DataFrame.from_records(H.prot_pep_pairs(), columns=('prot','pep'))
prot2pep['pep_cnt'] = prot2pep.pep.map(len)
prot_pep_counts     = prot2pep.groupby('prot').pep_cnt.sum()


def simple_postprocessing(prot_I, prot_reps, I_cols):
    prot_I = prot_I.join(prot_reps)
    cols = list(prot_reps.columns) + I_cols
    prot_I = prot_I[cols]
    prot_I = prot_I.set_index('repr')
    prot_I.index.name = 'protein'
    return prot_I

prot_deconvoluted = simple_postprocessing(prots_I, prot_reps, I_cols)
prot_deconvoluted.to_csv(Path(r"~/Projects/pep2prot/simple_protein_report.csv").expanduser())

prots_min_I['reps'] = prot_reps
prots_max_I['reps'] = prot_reps


