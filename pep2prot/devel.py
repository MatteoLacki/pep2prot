%load_ext autoreload
%autoreload 2

from collections import Counter
from networkx import Graph, connected_components
from pathlib import Path
import pandas as pd
from pandas import DataFrame
from functools import partial

pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', -1)#display whole column without truncation

path = Path(r"~/Projects/pep2prot/pep2prot/data/peptide_report.csv").expanduser()
D = pd.read_csv(path, encoding = "ISO-8859-1")

pep  = ['sequence', 'modifier'] # column defining a peptide
prot = ['accession'] # column names defining a protein
pep_prot = pep + prot

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


def map_edges_to_connected_components(S, pep_prot):
    """Map components to edges.

    Args:
        S (pandas.DataFrame): Run intensities indexed by pep-prot edges.
        pep_prot (list of str): Columns designating peptides and proteins.
    Returns:
        pandas.Series: maps pep-prots to connected components numbers.
    """
    G = Graph()
    G.add_edges_from(((s,m),a) for s,m,a in S.index)# (seq,mod,acc)
    edge2cc = ((a[0],a[1],b,i) if isinstance(a, tuple) else (b[0],b[1],a,i)
               for i, cc in enumerate(connected_components(G))
               for a, b in G.subgraph(cc).edges)
    return DataFrame(edge2cc, columns=pep_prot+['cc']).set_index(pep_prot).cc


S = aggregate_pep_clust(D, [c for c in D.columns if "intensity in" in c], pep_prot)
S['cc'] = map_edges_to_connected_components(S, pep_prot)
max_cc_idx = S.groupby('cc').size().argmax()

# pep_prot_edges_in_connected_components = Counter(S.groupby('cc').size())
S_cc = S.groupby('cc')
CC = S_cc.get_group(max_cc_idx)








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

