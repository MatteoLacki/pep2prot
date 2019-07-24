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
import re

from pep2prot.read import read_isoquant_peptide_report, read_fastas
from pep2prot.graphs import ProtPepGraph, BiGraph, get_peptide_protein_graph
from pep2prot.preprocessing import preprocess_isoquant_peptide_report, complex_cluster_buster, simple_cluster_buster, get_sequences
from pep2prot.intensities import get_prot_intensities
from pep2prot.postprocessing import get_info_on_prots, get_stats
from pep2prot.range_ops import range_list_len


min_pepNo_per_prot = 2
max_rt_deviation = 1
path = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()

D = read_isoquant_peptide_report(path/'peptide_report.csv')
D, I_cols = preprocess_isoquant_peptide_report(D)
unique_columns = ['peptide_overall_max_score','peptide_fdr_level','peptide_overall_replication_rate','prots','pre_homology_accessions','pi','mw']
# check that those do not change
# do I really need all of that???
observed_prots = {r for rg in D.prots for r in rg}
fastas = read_fastas(path/'mouse.fasta', observed_prots)


X = pd.DataFrame(((p,ps,s,e,r) 
                   for p,s,e,rg,ps in zip(D.index,
                                          D.start,
                                          D.end,
                                          D.prots,
                                          D.sequence) for r in rg),
                 columns=('pep','pepseq','start','end','prot'))
X['protseq'] = X.prot.map(fastas.prot_seq)
X = X.set_index(['pep','prot'])
assert all(ps in rs for ps,rs in zip(X.pepseq, X.protseq)), "Some peptides are not subsequences of proteins they are reported to explain."


X['protsubseq'] = [rs[(s-1):e] for s,e,rs in zip(X.start,X.end,X.protseq)]
X_ok  = X[X.protsubseq == X.pepseq].copy() # start-end point to the pepstring
X_bad = X[X.protsubseq != X.pepseq].copy() # on these, I HAVE to recalculate start-end
# on other, I MIGHT do it. And I should do it for comparison.

def starts_and_ends(subseq, seq):
    """Get the covered areas.

    fitting aaa to ABAaaaaB will result in [(3,6), (4,7)].
    """
    for m in re.finditer("(?=" + subseq + ")", seq):
        s = m.start()
        e = s + len(subseq)
        yield (s,e)

# long
correct_se = [frozenset(starts_and_ends(ps,rs)) for ps,rs in zip(X_ok.pepseq, X_ok.protseq)]
# short
rectified_se = [frozenset(starts_and_ends(*x)) for x in zip(X_bad.pepseq, X_bad.protseq)]


d = X_bad.groupby('prot').get_group(protseq)

for pepseq in d.pepseq:
    protseq = d.protseq[0]
    for se in starts_and_ends(pepseq, protseq):
        print(se)

s = []
for protseq,d in X_bad.groupby('prot'):
    s.append(sum_interval_lengths(se for pepseq in d.pepseq
                                for se in starts_and_ends(pepseq, d.protseq[0])))
starts_and_ends(*x) for x in zip(X_bad.pepseq, X_bad.protseq)

Counter(len(x) for x in correct_se)
Counter(len(x) for x in rectified_se)

X_ok['subseq_cnt'] = [len(x) for x in correct_se]
X_ok['found_se'] = correct_se
X_ok_ommited = X_ok.query('subseq_cnt != 1')
x = X_ok_ommited.iloc[0]

# need an example of prot: [(s,e),...,(s,e)]

def iter_ranges(X):
    for prot, d in X.groupby('prot'):
        protseq = d.protseq[0]
        yield prot, len(protseq), [se for pepseq in d.pepseq for se in starts_and_ends(pepseq, protseq) ]
# too many things: run drop_duplicates() somewhere
prot, prot_len, ranges = next(iter_ranges(X_bad))
x = {prot: range_list_len(ranges)/prot_len 
     for prot, prot_len, ranges in iter_ranges(X_bad)}
# numpyfying it all

X = pd.DataFrame(((p,ps,s,e,r) 
                   for p,s,e,rg,ps in zip(D.index,
                                          D.start,
                                          D.end,
                                          D.prots,
                                          D.sequence) for r in rg),
                 columns=('pep','pepseq','start','end','prot'))
X['protseq'] = X.prot.map(fastas.prot_seq)
X = X.drop_duplicates()

pepseqINprotseqCNT = [psp.count(p) for p, psp in zip(X.pepseq, X.protseq)]
Counter(pepseqINprotseqCNT)

i = next(i for i,j in enumerate(pepseqINprotseqCNT) if j == 10)
p, psp = X.iloc[i][['pepseq','protseq']]


from pep2prot.string_ops import find_indices3

%%time
W = pd.DataFrame.from_records((prot,find_indices3(pepseq,protseq)) 
                              for prot,protseq,pepseq in zip(X.prot,
                                                             X.protseq,
                                                             X.pepseq))


W = pd.DataFrame.from_records(((prot,s,e) for prot,protseq,pepseq in zip(X.prot,
                                                                         X.protseq,
                                                                         X.pepseq)
                                          for s,e in starts_and_ends(pepseq, protseq)),
                              columns=('prot','s','e'),
                              index='prot')
W = W.drop_duplicates()
W = W.sort_values(['prot', 's'])
prot_se0 = W.groupby('prot').head(1)

W.s - prot_se0.s
W.e - prot_se0.s

g,d = next((g,d) for g,d in iter(W_prot) if len(d)>10)

# d = pd.DataFrame({'s':[1,2,3,5,17,19,25], 'e':[10,9,8,11,19,100,27]})
s = np.array([1,2,3,5,17,19,25])
e = np.array([10,9,8,11,19,21,27])
# 1, 11 > 10
# 17,19 > 2
# 19,21 > 2
# 25,27 > 2
# ===========
#         16

s0,e0 = s[:-1],e[:-1]
s1,e1 = s[1:], e[1:]
cover_len = np.where(s1 >= e0, e1-s1, np.maximum(e1-e0,0))
cover_len.sum() + e0[0]-s0[0]
cover_len += 
# Fix, as it doesn't sum to what should be there







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

prots     = prots_I.index
prot_info = get_info_on_prots(prots, fastas)# find a better name for that function for God's sake...

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


