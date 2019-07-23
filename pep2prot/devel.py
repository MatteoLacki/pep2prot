%load_ext autoreload
%autoreload 2

from pathlib import Path
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import networkx as nx
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 30)#display whole column without truncation

from pep2prot.read import read_isoquant_peptide_report, read_fastas
from pep2prot.graphs import ProtPepGraph, BiGraph, get_peptide_protein_graph
from pep2prot.preprocessing import preprocess_isoquant_peptide_report, complex_cluster_buster, simple_cluster_buster, are_pepseqs_in_protseqs
from pep2prot.intensities import get_prot_intensities
from pep2prot.postprocessing import get_info_on_prots, get_stats

min_pepNo_per_prot = 2
max_rt_deviation = 1
path = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()

D = read_isoquant_peptide_report(path/'peptide_report.csv')
D, I_cols = preprocess_isoquant_peptide_report(D)
unique_columns = ['peptide_overall_max_score','peptide_fdr_level','peptide_overall_replication_rate','prots','pre_homology_accessions','pi','mw']
observed_prots = {r for rg in D.prots for r in rg}
fastas         = read_fastas(path/'mouse.fasta', observed_prots)
are_pepseqs_in_protseqs(D, fastas)
DD = complex_cluster_buster(D, I_cols, unique_columns, max_rt_deviation)
# DD = simple_cluster_buster(D, I_cols, unique_columns)
H, RWEP, BRR = get_peptide_protein_graph(DD) # aggregate the same peptides in different clusters
# H, RWEP, BRR = get_peptide_protein_graph(D3)
pep2pepgr = {p:pg for pg in H.peps() for p in pg}
DD        = DD.loc[pep2pepgr] # peps in H (no simple prot-pep pairs)
peps_I    = DD[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities
prots_min_I, prots_I, prots_max_I = get_prot_intensities(H, peps_I)
# prots_stats = get_stats(prots_min_I, prots_I, prots_max_I)

prots     = prots_I.index
prot_info = get_info_on_prots(prots, fastas)# find a better name for that function for God's sake...

prot2pep            = pd.DataFrame.from_records(H.prot_pep_pairs(), columns=('prot','pep'))
prot2pep['pep_cnt'] = prot2pep.pep.map(len)
prot_pep_counts     = prot2pep.groupby('prot').pep_cnt.sum()
# calculating the coverage!


prot2pepLong = pd.DataFrame(((rg, p) for rg, pg in H.prot_pep_pairs() 
                                     for p  in pg),
                            columns=('prot','pep')).set_index('prot')

prot2pepLong = prot2pepLong.join(D[['start','end']], on='pep')


s0,e0 = strano_start.loc['AGVHIK'][['start','end']].iloc[1,:]
s1,e1 = strano_start.loc['AGVHIK'][['start','end']].iloc[0,:]
KCRU_MOUSE_seq = fastas.loc['KCRU_MOUSE'].sequence
KCRB_MOUSE_seq = fastas.loc['KCRB_MOUSE'].sequence

KCRU_MOUSE_seq[s0-1:e0]
KCRU_MOUSE_seq[s1-1:e1]
KCRB_MOUSE_seq[s1-1:e1]



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

# OK, now need to add in some protein specific information?


