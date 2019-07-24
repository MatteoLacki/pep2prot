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
max_rt_deviation = 1
path = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()

D = read_isoquant_peptide_report(path/'peptide_report.csv')
D, I_cols = preprocess_isoquant_peptide_report(D)
unique_columns = ['peptide_overall_max_score','peptide_fdr_level','peptide_overall_replication_rate','prots','pre_homology_accessions','pi','mw']

fastas = read_fastas(path/'mouse.fasta', {r for rg in D.prots for r in rg})
prots  = get_protein_coverages(D, fastas)
# plt.hist(prots.pep_coverage, bins=100)
# plt.show()

DD = complex_cluster_buster(D, I_cols, unique_columns, max_rt_deviation)
# DD = simple_cluster_buster(D, I_cols, unique_columns)
H, RWEP, BRR = get_peptide_protein_graph(DD) # aggregate the same peptides in different clusters
# H, RWEP, BRR = get_peptide_protein_graph(D3)

pep2pepgr = {p:pg for pg in H.peps() for p in pg}
DDinH     = DD.loc[pep2pepgr] # peps in H: no simple prot-pep pairs, no unnecessary prots?
DDinH['pepgr'] = DDinH.index.map(pep2pepgr)
peps_I    = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities

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
# plt.hist(prots_HF.pep_coverage, bins=100)
# plt.show()



prot_info = get_info_on_prots(prots_I.index, fastas)# find a better name for that function for God's sake...

prot2pep            = pd.DataFrame.from_records(H.prot_pep_pairs(), columns=('prot','pep'))
prot2pep['pep_cnt'] = prot2pep.pep.map(len)
prot_pep_counts     = prot2pep.groupby('prot').pep_cnt.sum()
# calculating the coverage!
# restricting to observed ones
pepseq2protseq = pepseq2protseq.loc[prot_info['representative protein']]


# 42553 peptides are assigned to these bloody proteins.
len([pg for pg in H.peps() for p in pg])

T = pepseq2protseq.iloc[0:100]
pepseq, protseq = T.iloc[0]



x = [list(iter_starts_and_ends(pepseq, protseq)) for pepseq, protseq in zip(T.pepseq.values, T.protseq.values)]
# this is fast for small cases.

# IT should be like that: if I want to caclulate this, ok.
# Also, if the suggested position doesn't match, calculate it.
# Or, if you don't believe PLGS, calculate everything.

# So, unfortunately I still need the boolean function that checks, if the suggested positions are correct.

# 0. get the starts and ends of peptides in H

# differently: there are mostly how many prots per pep??
D['prots_cnt'] = D.prots.map(len)
start_end_unique = D[D.prots_cnt == 1]
start_end_unique.sequence
start_end_unique.join(fastas, on='prots_str')


start_end_nonunique = D[D.prots_cnt != 1]




protsPERpep = Counter() # mostly one protein per peptide
sum(protsPERpep.values()) - protsPERpep[1]
# the problem is, that start and end might correspond to many things if there are more proteins.
D[['start','end']]



    



prot2pepLong = pd.DataFrame(((rg, p) for rg, pg in H.prot_pep_pairs() 
                                     for p  in pg),
                            columns=('prot','pep')).set_index('prot')




# prot2pepLong = prot2pepLong.join(D[['start','end']], on='pep')

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


