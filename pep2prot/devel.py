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
from pep2prot.postprocessing import summarize_prots, get_stats, prettify_protein_informations, get_full_report

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
pep2pepgr      = {p:pg for pg in H.peps() for p in pg}
DDinH          = DD.loc[pep2pepgr] # peps in H: no simple prot-pep pairs, no unnecessary prots?
DDinH['pepgr'] = DDinH.index.map(pep2pepgr)
peps_I         = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities
prots_min_I, prots_I, prots_max_I = get_prot_intensities(H, peps_I)
prot_info      = summarize_prots(H, fastas, prots.pep_coverage)
prots_I_nice   = prettify_protein_informations(prots_I, prot_info)

prots_I_nice.to_csv(Path(r"~/Projects/pep2prot/protein_report.csv").expanduser())
# prots_stats = get_stats(prots_min_I, prots_I, prots_max_I)

all_prots = get_full_report(prots_min_I, prots_I, prots_max_I)
add_info_to_intensities_nice = prettify_protein_informations(all_prots, prot_info)
add_info_to_intensities_nice.to_csv(Path(r"~/Projects/pep2prot/full_protein_report.csv").expanduser())

# DEBUG
# prots_I
# error_maker = 'MPRIP_MOUSE'
# sample = "intensity in 2019-023-22 Wolf 1"
# ergr = [rg for rg in prots_I.index for r in rg if error_maker in r][0]
# er_peps = [p for pg in H[ergr] for p in pg]

# 35534/2
# H[frozenset(['IEDLQR'])]
# H.subgraph(nx.node_connected_component(H, ergr))

# def peep():
#     for pg in H[ergr]:
#         yield pg
#         for rg in H[pg]:
#             yield rg
# MARKOV_BLANKET = H.subgraph(set(peep()))
# MARKOV_BLANKET.draw(with_labels=True)

# D.loc[er_peps,sample]
# the_other_protein = list(MARKOV_BLANKET.prots())[1]
# prots_I_nice.loc[error_maker,sample]
# prots_I_nice.loc[the_other_protein,sample]