%load_ext autoreload
%autoreload 2

import matplotlib.pyplot as plt
from collections import Counter
from networkx import Graph, connected_components
import networkx as nx
import numpy as np
from pathlib import Path
import pandas as pd
from pandas import DataFrame

pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', -1)#display whole column without truncation

path = Path(r"~/Projects/pep2prot/pep2prot/data/peptide_report.csv").expanduser()

D = pd.read_csv(path, encoding = "ISO-8859-1")
# preprocessing
D.modifier = D.modifier.fillna('')
I_cols = [c for c in D.columns if "intensity in" in c]
D[I_cols] = D[I_cols].fillna(0)
D.rename(columns={'pre_homology_entries':'prots'}, inplace=True)
D['pep'] = np.where(D['modifier'], D['sequence'] + "_" + D['modifier'], D['sequence'])
pep_prots = ['pep','prots']
# semplifying
pep = D.pep
prots = D.prots
G = Graph()# peptides start with a $
G.add_edges_from((("$"+pe,p) for pe, pr in zip(pep, prots) for p in pr.split(',')))

def nodes_poor_in_neighbors(G, k, node_selector):
    for n in G:
        if node_selector(n) and G.degree(n) < k:
            yield n
is_peptide = lambda p: p[0]=='$'
is_protein = lambda p: p[0]!='$'

def count_peptides_proteins(G):
    pe = pr = 0
    for n in G:
        if is_protein(n):
            pr += 1
        else:
            pe += 1
    return pr, pe

print(count_peptides_proteins(G))
G.remove_nodes_from(list(nodes_poor_in_neighbors(G, 2, is_protein)))
print(count_peptides_proteins(G))
G.remove_nodes_from(list(nodes_poor_in_neighbors(G, 1, is_peptide)))
print(count_peptides_proteins(G))

Counter(count_peptides_proteins(G.subgraph(cc)) for cc in nx.connected_components(G))



S = D.groupby('pep')[I_cols].sum()
S['prots'] = [",".join(G["$"+p]) for p in S.index]



edge2cc = ((a[1:],b,c) if a[0] == '$' else (b[1:],a,c)
           for c,cc in enumerate(connected_components(G))
           for a, b in G.subgraph(cc).edges)# c = number of con-comp cc
cc = DataFrame(edge2cc, columns=['pep','prot','cc'])
CC = list(connected_components(G))







cc_summary = cc.groupby('cc').nunique()
cc.query('cc == 4')


cc_summary.groupby(['pep','prot']).size()
W = cc_summary.groupby(['prot','pep']).size().reset_index()
dict(W.groupby('prot').size()).items()

# the big cluster???


# .set_index(pep_prot).cc
# if return_graph:
#     return cc, G
# else:
#     return cc





#TODO: add in some filtering here
def aggregate_pep_clust(D, intensity_columns, pep_prot):
    """Sum up peptide intensities found across different clusters.

    Attention: if in a given run a peptide was not found, we assign it zero intensity.

    Args:
        D (pandas.DataFrame): IsoQuant report from csv.
        intensity_columns (list of str): Columns containing run intensities in D.
        pep_prot (list of str): Columns designating peptides and proteins.
    """
    return D.groupby(pep_prot)[intensity_columns].sum()


def map_edges_to_connected_components(S, pep_prot, return_graph=False):
    """Map components to edges.

    Args:
        S (pandas.DataFrame): Run intensities indexed by pep-prot edges.
        pep_prot (list of str): Columns designating peptides and proteins.
        return_graph (boolean): Return the intermediate pep-prot graph?
    Returns:
        pandas.Series: maps pep-prots to connected components numbers.
    """
    G = Graph()
    G.add_edges_from(((s,m),a) for s,m,a in S.index)# (seq,mod,acc)
    edge2cc = ((a[0],a[1],b,c) if isinstance(a, tuple) else (b[0],b[1],a,c)
               for c,cc in enumerate(connected_components(G))
               for a, b in G.subgraph(cc).edges)# c = number of con-comp cc
    cc = DataFrame(edge2cc, columns=pep_prot+['cc']).set_index(pep_prot).cc
    if return_graph:
        return cc, G
    else:
        return cc

S = aggregate_pep_clust(D, [c for c in D.columns if "intensity in" in c], pep_prot)
S['cc'], G = map_edges_to_connected_components(S, pep_prot, True)
max_cc_idx = np.argmax(S.groupby('cc').size().values)


def count_proteins_per_ccs(G):
    PP = Counter()
    for i, cc in enumerate(connected_components(G)):
        for n in cc:
            if not isinstance(n, tuple):
                PP[i] += 1
    return Counter(PP.values())

from itertools import islice
import matplotlib.pyplot as plt

count_proteins_per_ccs(G)

k = 100
N_k = set.union(*islice(connected_components(G), k))
nx.draw(G.subgraph(N_k))
plt.show()


# pep_prot_edges_in_connected_components = Counter(S.groupby('cc').size())
S_cc = S.groupby('cc')
CC = S_cc.get_group(max_cc_idx)

# count the peps and prots in different graphs
W.groupby('cc').get_group(max_cc_idx).accession.unique()
W = S.reset_index()
W['pep'] = W['sequence'] + "_" + W['modifier']
peps_prots_per_cc = W.groupby('cc')[prot+['pep']].nunique()
peps_prots_summary = peps_prots_per_cc.groupby(prot+['pep']).size()
# do I do something wrong, or there are no multiple proteins per cc?







# def get_peptide_potein_graph(report, intensity_estimator,
#                              group_vars=['sequence', 'modifier'],
#                              prot_tag="accession"):
#     """Get the peptide protein graph."""
#     G = nx.Graph()
#     for peptide, pep_df in report.groupby(group_vars):
#         G.add_node(peptide, I=intensity_estimator(pep_df))
#         for protein in set(pep_df[prot_tag]):
#             G.add_edge(peptide, protein)
#     return G

# # def total_signal_intensity(peptide_df):
# #     return peptide_df.signal_intensity.sum()

# # def top3_signal_intensity(peptide_df):
# #     return peptide_df.signal_intensity.nlargest(3).sum()

# def total_intensity(peptide_df):
    
#     return peptide_df[I_cols].sum()

# G = get_peptide_potein_graph(report, total_intensity)
# len(G)


# CCS = list(nx.connected_components(G))
# CCS_counts_of_node_numbers = Counter(len(cc) for cc in CCS)

# Counter(report.groupby(['sequence', 'modifier']).size())





# def deconvolve_intensities(connected_pep2prot_graph):
#     pass

