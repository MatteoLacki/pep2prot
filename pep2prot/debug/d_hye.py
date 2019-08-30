%load_ext autoreload
%autoreload 2

import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 30)#display whole column without truncation
from pandas import DataFrame as df
from pathlib import Path
from collections import Counter, defaultdict
import networkx as nx
import matplotlib.pyplot as plt
from plotnine import *

from pep2prot.read import read_isoquant_peptide_report, read_fastas
from pep2prot.preprocessing import preprocess_isoquant_peptide_report, get_protein_coverages, complex_cluster_buster, simple_cluster_buster 
from pep2prot.graphs.pep_prot_graph import ProtPepGraph
from pep2prot.intensities import get_prot_intensities
from pep2prot.df_ops import sum_top
from pep2prot.postprocessing import summarize_prots, get_stats, prettify_protein_informations, get_full_report


test_data = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()
pep_rep_path = test_data/'hye_peprep.csv'
fastas_path = test_data/'HYE.fasta'
cluster_buster = complex_cluster_buster

D = read_isoquant_peptide_report(pep_rep_path)
D, I_cols = preprocess_isoquant_peptide_report(D)

fastas = read_fastas(fastas_path)
observed_prots = {r for rg in D.prots for r in rg}
assert all(r in fastas.index for r in observed_prots), "It seems that you are using a different set of fastas than the peptide annotation software before. Repent please."
prots = get_protein_coverages(D, fastas)
uni_cols = ['peptide_overall_max_score','peptide_fdr_level',
        'peptide_overall_replication_rate','prots',
            'pre_homology_accessions','pi','mw']

DD = cluster_buster(D, I_cols, uni_cols) # agg same peptides in various clusters
assert DD.index.is_unique, "Some peptides are still across different clusters."

G = ProtPepGraph((r,p) for p, rg in zip(DD.index, DD.prots) for r in rg)
lonely, unsupported = G.remove_lonely_and_unsupported(2)
H, rejected = G.get_minimal_graph()

pep2pepgr = {p:pg for pg in H.peps() for p in pg}
DDinH = DD.loc[pep2pepgr] # peps in H: no simple prot-pep pairs, no unnecessary prots?
DDinH['pepgr'] = DDinH.index.map(pep2pepgr)

# implementing the max three peptide rule: accross runs and inside runs
pep_I = DDinH[I_cols]

# (prot_cnt, pep_cnt), pep2prot_cnt = H.nodes_cnt(), len(H.edges)
_, Icols_cnt = pep_I.shape
Icols = pep_I.columns

pep2prot = pd.MultiIndex.from_tuples((R,p) for R,P in H.prot_pep_pairs() for p in P)
pep2prot.names = ('prot','pep')
pep2prot_maxI = df(index=pep2prot).join(pep_I, on='pep')

prot_maxI = sum_top(pep2prot_maxI.droplevel('pep'), 3, 'prot')

unipep = set(p for P in H.peps(deg=1) for p in P)
unipep2prot_I = pep2prot_maxI[pep2prot.isin(unipep, 'pep')]
uniprot_minI = sum_top(unipep2prot_I.droplevel('pep'), 3, 'prot')

uniprot = uniprot_minI.index
inuprot = {R for R in H.prots() if not R in uniprot}
inuprot_noI = df(np.zeros(shape=(len(inuprot), Icols_cnt)),
                          index=inuprot,
                          columns=Icols)




# prot = uniprot ⊔ inuprot
prot_minI = pd.concat([uniprot_minI, inuprot_noI])

def get_meds(X):
    meds = pd.concat([X.filter(regex='HYE Mix A').median(axis=1),
                      X.filter(regex='HYE Mix B').median(axis=1)],
                      axis=1)
    meds.columns = ['A', 'B']
    return meds

def get_means(X):
    meds = pd.concat([X.filter(regex='HYE Mix A').mean(axis=1),
                      X.filter(regex='HYE Mix B').mean(axis=1)],
                      axis=1)
    meds.columns = ['A', 'B']
    return meds

pep_I_old = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities
prot_minI_old, prot_I_old, prot_maxI_old = get_prot_intensities(H, pep_I_old)

# X = get_meds(prot_minI)
# Y = get_meds(prot_minI_old)
X = get_means(prot_minI)
Y = get_means(prot_minI_old)

X['algo'] = 'top3'
Y['algo'] = 'all intensities'
Z = pd.concat([X,Y])

(ggplot(Z) + 
 geom_hline(yintercept=np.log([1/2,1,4]))+
 geom_point(aes(x='np.log(A)', y='np.log(B)-np.log(A)'), size=2) + 
 facet_grid('.~algo', scales='free_x') )

# investigating sample-free intensity convergences
X = prot_minI_old.filter(regex='HYE Mix A').copy()
X.columns = ['a','b','c']
X = pd.concat([X.a, np.log(X.b) - np.log(X.a), np.log(X.c) - np.log(X.a)], axis=1)
X = X[X.a != 0]
X.columns = ['a', 'b2a', 'c2a']
plt.scatter(np.log(X.a), X.b2a, s=1)
plt.scatter(np.log(X.a), X.c2a, s=1)
plt.show()


X = prot_minI_old.filter(regex='HYE Mix A').copy()
X.columns = ['a','b','c']

plt.scatter(x=np.log(X.a), y=X.std(axis=1), s=1)
plt.show()

plt.scatter(x=np.log(X.median(axis=1)), y=X.mad(axis=1), s=1)
plt.show()

(ggplot(Y) + 
 geom_point(aes(x='np.log(A)', y='np.log(B)-np.log(A)'), size=.2) )

plt.scatter(x=np.log(X.median(axis=1)), y=np.log(X.mad(axis=1)), s=1)
plt.show()

# comparing the values for plots of peptide I, rather than protein I
A = prot_I_old.filter(regex='Mix A')
B = prot_I_old.filter(regex='Mix B')
A_med = A.median(axis=1)
B_med = B.median(axis=1)

p = np.linspace(50,100, 100)
q = np.percentile(A.median(axis=1), q=p)
# plt.plot(q,p)
# plt.show()
q = np.percentile(A_med, q=95)

HI_prots = {r for R in A_med[A_med >= q].index for r in R}
HI_prots_species = defaultdict(set)
for r in HI_prots:
    acc, species = r.split('_')
    HI_prots_species[species].add(species)

 
HI_peps = {p for p, R in zip(DD.index, DD.prots) if any(r in HI_prots for r in R)}
HI_peps_I = DD.loc[HI_peps,I_cols]

Z = get_means(HI_peps_I)
plt.scatter(np.log(Z.A), np.log(Z.A)-np.log(Z.B))
plt.show()
