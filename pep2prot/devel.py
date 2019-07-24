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
DD     = cluster_buster(D, I_cols, unique_columns)# agg same peptides in various clusters
assert np.all(DD.groupby(['pep','prots']).size() == 1), "Some peptides are still across different clusters."
H, prots_no_peps, peps_no_prots, beckham_prots = get_peptide_protein_graph(DD) 
pep2pepgr = {p:pg for pg in H.peps() for p in pg}
DDinH           = DD.loc[pep2pepgr] # peps in H: no simple prot-pep pairs, no unnecessary prots?
DDinH['pepgr']  = DDinH.index.map(pep2pepgr)
peps_I          = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities
prots_min_I, prots_I, prots_max_I = get_prot_intensities(H, peps_I)
# prots_stats = get_stats(prots_min_I, prots_I, prots_max_I)
# plt.hist(prots_HF.pep_coverage, bins=100); plt.show()
prot_info = get_info_on_prots(prots_I.index, fastas)# find a better name for that function for God's sake...
prot_info = prot_info.join(prots.pep_coverage, on='representative protein')
prot2pep = pd.DataFrame.from_records(H.prot_pep_pairs(), columns=('prot','pep'))
prot2pep['pep_cnt'] = prot2pep.pep.map(len)
prot_info['peptides count'] = prot2pep.groupby('prot').pep_cnt.sum()




def add_info_to_intensities(prot_info, prot_intensities):
    X = prot_info.join(prot_intensities.applymap(int))
    X = X.set_index('representative protein')
    X.pep_coverage = ["{:.3%}".format(pc) for pc in X.pep_coverage]
    X = X.rename(columns={'prot_seq':'sequence', 'seq_len':'sequence length', 'pep_cover': 'peptide coverage', 'protein_group': 'other proteins in group'})
    X = X.drop(columns='sequence')
    X['monoisotopic mass'] = [round(m,3) for m in X['monoisotopic mass']]
    return X

prots_I_nice = add_info_to_intensities(prot_info, prots_I)
prots_I_nice.to_csv(Path(r"~/Projects/pep2prot/protein_report.csv").expanduser())

all_prots = np.zeros(dtype=prots_I.values.dtype, shape=(prots_I.shape[0],prots_I.shape[1]*3))
all_prots[:,0::3] = prots_min_I.values
all_prots[:,1::3] = prots_I.values
all_prots[:,2::3] = prots_max_I.values

columns = []
for I in prots_min_I.columns:
    columns.append('minimal '+I)
    columns.append('deconvoled '+I)
    columns.append('maximal '+I)    
all_prots = pd.DataFrame(all_prots, columns=columns, index=prots_I.index)
add_info_to_intensities_nice = add_info_to_intensities(prot_info, all_prots)
add_info_to_intensities_nice.to_csv(Path(r"~/Projects/pep2prot/full_protein_report.csv").expanduser())