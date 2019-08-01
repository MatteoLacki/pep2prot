%load_ext autoreload
%autoreload 2

import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 30)#display whole column without truncation
from pathlib import Path
from collections import Counter
import networkx as nx

from pep2prot.read import read_isoquant_peptide_report, read_fastas
from pep2prot.graphs import ProtPepGraph, get_full_prot_pep_graph, BiGraph
from pep2prot.preprocessing import preprocess_isoquant_peptide_report, get_protein_coverages, complex_cluster_buster, simple_cluster_buster 
from pep2prot.intensities import get_prot_intensities
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
assert np.all(DD.groupby(['pep','prots']).size() == 1), "Some peptides are still across different clusters."

G = ProtPepGraph((r,p) for p, rg in zip(DD.index, DD.prots) for r in rg)
lonely, unsupported = G.remove_lonely_and_unsupported(2)
H, beckham_prot_groups = G.get_minimal_graph()

pep2pepgr = {p:pg for pg in H.peps() for p in pg}
DDinH = DD.loc[pep2pepgr] # peps in H: no simple prot-pep pairs, no unnecessary prots?
DDinH['pepgr'] = DDinH.index.map(pep2pepgr)
peps_I = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities
peps_I.isna().any().any()
debug = True
# go to intensities.py


buggy_weights = weights[weights.isnull().any(axis=1)]
buggy_pep_groups = set(buggy_weights.index.get_level_values('pep'))
buggy_peps = {p for pg in buggy_pep_groups for p in pg}
D.loc[buggy_peps]
buggy_prot_groups = {rg for pg in buggy_pep_groups for rg in H[pg]}
# H.subgraph(buggy_pep_groups|buggy_prot_groups).draw(with_labels=True) # looks normal

rg = next(iter(set(buggy_weights.index.get_level_values('prot'))))
H.subgraph(nx.node_connected_component(H, rg)).draw(with_labels=True, font_size=6)




weights.xs(pep=pg, prot=rg)

weights.loc[pg]
weights.index
pg in weights.index.get_level_values('pep')
rg in weights.index.get_level_values('prot')
(pg, rg) in weights.index

prots_min_I, prots_I, prots_max_I = get_prot_intensities(H, peps_I)

weights.loc[buggy_pep_groups]
x = weights.xs(pg, level='pep')
x = x.isna().any(axis=1)
rg = next(iter(x[x].index))
D.loc['LAADDFR', I_cols]
# there might be a need for left or right merge somewhere

weights.xs((pg,rg), level=('pep','prot'))
weights
prots_curr_I.loc[rg]
prots_curr_I.loc[]
prots_curr_I.loc[[frozenset({'K1C15_HUMAN', 'K1C15_HUMAN_CONTA'})]]
prots_curr_I.loc[[rg]]
# OK, so there is a bug! some protein group is simply not there. how come???
set(otherpeps2mixprots.get_level_values('prot')) - set(prots_curr_I.index)
H[rg]
set(peps2prots.get_level_values('pep')) - set(otherpeps_I.index)

# in principle, it works.
edges = [('A',0), ('A',1), ('B',1), ('B',2), ('B',4), ('C',2), ('C',3), ('C',4), ('C',5), ('D',4), ('D',5), ('D',6)]
test = ProtPepGraph(edges)
H_test, beckham_prot_groups_test = test.get_minimal_graph()
H_test.draw(with_labels=True)


get_full_prot_pep_graph()


prot_info      = summarize_prots(H, fastas, prots.pep_coverage)
prots_I_nice   = prettify_protein_informations(prots_I, prot_info)
prots_stats    = get_stats(prots_min_I, prots_I, prots_max_I)

if verbose:
    print('Preparing reports.')
all_prots      = get_full_report(prots_min_I, prots_I, prots_max_I)
all_prots_nice = prettify_protein_informations(all_prots, prot_info)




