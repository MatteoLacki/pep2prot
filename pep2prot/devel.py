%load_ext autoreload
%autoreload 2

import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from networkx import Graph, connected_components
import networkx as nx
import numpy as np
from pathlib import Path
import pandas as pd
from pandas import DataFrame
from itertools import islice

from pep2prot.pep_prot_graph import PepProtGraph
from pep2prot.bipartite_graph import BiGraph

pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', -1)#display whole column without truncation
path = r"~/Projects/pep2prot/pep2prot/data/peptide_report.csv"
path = Path(path).expanduser()
D = pd.read_csv(path, encoding = "ISO-8859-1")
verbose = True
peps_per_prot = 2

# preprocessing
D.modifier = D.modifier.fillna('')
I_cols = [c for c in D.columns if "intensity in" in c]
D[I_cols] = D[I_cols].fillna(0)
D.rename(columns={'pre_homology_entries':'prots'}, inplace=True)
D['pep'] = np.where(D['modifier'], D['sequence'] + "_" + D['modifier'], D['sequence'])
pep_prots = ['pep','prots']
D.prots = D.prots.str.split(',')

# G = PepProtGraph(zip(D.pep, D.prots))
# G.del_prots_with_little_peps(2)
# prothoods = defaultdict(list)
# for p in G._peps(): #can there be peptide duplicates? No.
#     prot_neighbors = frozenset(G[p])
#     prothoods[prot_neighbors].append(p[1:])
# F = BiGraph((r, " | ".join(ps)) for rs, ps in prothoods.items() for r in rs)
# # again!


# # simplifying the graph
# G = BiGraph((r,p) for rs, p in zip(D.prots, D.pep) for r in rs)
# prothoods = defaultdict(list)
# for p in G.B():
#     prothoods[frozenset(G[p])].append(p)
# F = BiGraph((r, frozenset(ps)) for rs, ps in prothoods.items() for r in rs)
# pepgrouphoods = defaultdict(list)
# for r in F.A():
#     pepgrouphoods[frozenset(F[r])].append(r)
# H = BiGraph((frozenset(rs), ps) for ps, rs in pepgrouphoods.items())

# print(len(H))
# print(len(list(H.A())))
# print(len(list(H.B())))

def merge_nodes(G, nodes2merge):
    neighbors_of_nodes_to_merge = defaultdict(set)
    for n in nodes2merge:
        neighbors_of_nodes_to_merge[frozenset(G[n])].add(n)
    return BiGraph((m, frozenset(merged))
                   for to_merge, merged in neighbors_of_nodes_to_merge.items()
                   for m in to_merge)

G = BiGraph((r,p) for rs, p in zip(D.prots, D.pep) for r in rs)
# filter out proteins with not enough peptides
min_pepNo_per_prot = 2
prots_deg_01 = [r for r in G.A() if G.degree(r) < min_pepNo_per_prot]
G.remove_nodes_from(prots_deg_01)
peps_for_prots_deg_01 = [p for p in G.B() if G.degree(p) == 0]
G.remove_nodes_from(peps_for_prots_deg_01)

peps = G.B()
F = merge_nodes(G, peps)
prots = F.A()
H = merge_nodes(F, prots)

next(H.A()) # pep groups
next(H.B()) # prot groups
len(H) # graph is much smaller now
len(list(nx.connected_components(H)))

CC = next(nx.connected_components(H))
list(H.subgraph(CC).A())
list(H.subgraph(CC).B())



good = []
bad = []
for CC in nx.connected_components(H):
    CC = H.subgraph(CC)
    pepNo, protNo = CC.nodes_cnt()
    if protNo == 1:
        good.append(CC)
    else:
        bad.append(CC)

Counter(cc.nodes_cnt() for cc in good)# we can already iterate over these out to the output
Counter(cc.nodes_cnt() for cc in bad)# these need further analysis

cc = next(iter(bad)).copy() # this will not work otherwise?
# for cc in bad:
cc.nodes()
# cc.draw()

supported_prot_groups = [r for p in cc.A() if cc.degree(p) == 1 for r in cc[p]]

cc.remove_nodes_from(supported_prot_groups)
# just make the operation global: delete all those peptides.
# no need to work on the connected components now.



# I need the input fasta file to retrieve their sequence to choose the representative.


# it = nx.connected_components(H)
# H.subgraph(next(it)).draw()

# # iterating peps with deg=1 and their proteins
# good = list( (p,r) for p in G._peps() if G.degree(p) == 1 for r in G[p])
# good_peps, good_prots = map(set, list(zip(*good)))

# G.remove_nodes_from(good_peps)
# G.pep_degs()

# G.remove_nodes_from(good_prots)
# G.pep_degs()

# # prots that don't have anything now:
# G.remove_nodes_from(list(G.low_deg_nodes(0, G.is_prot)))

# # remove peps that do not give any new info:
# G.remove_nodes_from(list(G.low_deg_nodes(0, G.is_pep)))

# len(list(nx.connected_components(G)))

# it = nx.connected_components(G)
# n = next(it)
# W = G.subgraph(n)
# print(len(list(W.prots())))
# W.draw()

# G





# # This is mostly for reporting. Work on graphs first.
# S = D.groupby('pep')[I_cols].sum()
# S['prots'] = 

# [",".join(G["$"+p]) for p in S.index]

# edge2cc = ((a[1:],b,c) if a[0] == '$' else (b[1:],a,c)
#            for c,cc in enumerate(connected_components(G))
#            for a, b in G.subgraph(cc).edges)# c = number of con-comp cc
# cc = DataFrame(edge2cc, columns=['pep','prot','cc'])
# CC = list(connected_components(G))







# cc_summary = cc.groupby('cc').nunique()
# cc.query('cc == 4')


# cc_summary.groupby(['pep','prot']).size()
# W = cc_summary.groupby(['prot','pep']).size().reset_index()
# dict(W.groupby('prot').size()).items()

# # the big cluster???


# # .set_index(pep_prot).cc
# # if return_graph:
# #     return cc, G
# # else:
# #     return cc


# #TODO: add in some filtering here
# def aggregate_pep_clust(D, intensity_columns, pep_prot):
#     """Sum up peptide intensities found across different clusters.

#     Attention: if in a given run a peptide was not found, we assign it zero intensity.

#     Args:
#         D (pandas.DataFrame): IsoQuant report from csv.
#         intensity_columns (list of str): Columns containing run intensities in D.
#         pep_prot (list of str): Columns designating peptides and proteins.
#     """
#     return D.groupby(pep_prot)[intensity_columns].sum()


# def map_edges_to_connected_components(S, pep_prot, return_graph=False):
#     """Map components to edges.

#     Args:
#         S (pandas.DataFrame): Run intensities indexed by pep-prot edges.
#         pep_prot (list of str): Columns designating peptides and proteins.
#         return_graph (boolean): Return the intermediate pep-prot graph?
#     Returns:
#         pandas.Series: maps pep-prots to connected components numbers.
#     """
#     G = Graph()
#     G.add_edges_from(((s,m),a) for s,m,a in S.index)# (seq,mod,acc)
#     edge2cc = ((a[0],a[1],b,c) if isinstance(a, tuple) else (b[0],b[1],a,c)
#                for c,cc in enumerate(connected_components(G))
#                for a, b in G.subgraph(cc).edges)# c = number of con-comp cc
#     cc = DataFrame(edge2cc, columns=pep_prot+['cc']).set_index(pep_prot).cc
#     if return_graph:
#         return cc, G
#     else:
#         return cc

# S = aggregate_pep_clust(D, [c for c in D.columns if "intensity in" in c], pep_prot)
# S['cc'], G = map_edges_to_connected_components(S, pep_prot, True)
# max_cc_idx = np.argmax(S.groupby('cc').size().values)


# def count_proteins_per_ccs(G):
#     PP = Counter()
#     for i, cc in enumerate(connected_components(G)):
#         for n in cc:
#             if not isinstance(n, tuple):
#                 PP[i] += 1
#     return Counter(PP.values())

# from itertools import islice
# import matplotlib.pyplot as plt

# count_proteins_per_ccs(G)

# k = 100
# N_k = set.union(*islice(connected_components(G), k))
# nx.draw(G.subgraph(N_k))
# plt.show()


# # pep_prot_edges_in_connected_components = Counter(S.groupby('cc').size())
# S_cc = S.groupby('cc')
# CC = S_cc.get_group(max_cc_idx)

# # count the peps and prots in different graphs
# W.groupby('cc').get_group(max_cc_idx).accession.unique()
# W = S.reset_index()
# W['pep'] = W['sequence'] + "_" + W['modifier']
# peps_prots_per_cc = W.groupby('cc')[prot+['pep']].nunique()
# peps_prots_summary = peps_prots_per_cc.groupby(prot+['pep']).size()
# # do I do something wrong, or there are no multiple proteins per cc?







# # def get_peptide_potein_graph(report, intensity_estimator,
# #                              group_vars=['sequence', 'modifier'],
# #                              prot_tag="accession"):
# #     """Get the peptide protein graph."""
# #     G = nx.Graph()
# #     for peptide, pep_df in report.groupby(group_vars):
# #         G.add_node(peptide, I=intensity_estimator(pep_df))
# #         for protein in set(pep_df[prot_tag]):
# #             G.add_edge(peptide, protein)
# #     return G

# # # def total_signal_intensity(peptide_df):
# # #     return peptide_df.signal_intensity.sum()

# # # def top3_signal_intensity(peptide_df):
# # #     return peptide_df.signal_intensity.nlargest(3).sum()

# # def total_intensity(peptide_df):
    
# #     return peptide_df[I_cols].sum()

# # G = get_peptide_potein_graph(report, total_intensity)
# # len(G)


# # CCS = list(nx.connected_components(G))
# # CCS_counts_of_node_numbers = Counter(len(cc) for cc in CCS)

# # Counter(report.groupby(['sequence', 'modifier']).size())





# # def deconvolve_intensities(connected_pep2prot_graph):
# #     pass

