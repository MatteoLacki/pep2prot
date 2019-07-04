%load_ext autoreload
%autoreload 2

from collections import Counter
import networkx as nx
from pathlib import Path
import pandas as pd

pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', -1)#display whole column without truncation

path = Path(r"~/Projects/pep2prot/data/peptide_report.csv").expanduser()

def read_report(path):
    return pd.read_csv(path, encoding = "ISO-8859-1")

report = read_report(path)
I_cols = [c for c in report.columns if "intensity in" in c]

# p = report.groupby(['sequence', 'modifier', 'signal_z']).size()
# p[p>1]
# 989 / 8231 strange cases
# peptide_df = report.query("sequence == 'AAAHHYGAQCDK'")

# seq, mod = "AAAHHYGAQCDK", "Carbamidomethyl C(10)"
# d = report_pep_groups.get_group((seq, mod))

# First simplify with C++ code, then build graph!


simpler_report = report.groupby(['sequence', 'modifier', 'accession'])[I_cols].sum()

def pep_prot_iter(simpler_report):
    for seq, mod, acc in simpler_report.index:
        yield (seq, mod), acc

GG = nx.Graph()
GG.add_edges_from(pep_prot_iter(simpler_report))

def iter_edge2ccNo(GG):
    for i, cc in enumerate(nx.connected_components(GG)):
        for a, b in GG.subgraph(cc).edges:
            yield (a[0],a[1],b,i) if isinstance(a, tuple) else (b[0],b[1],a,i)
             
CCS = pd.DataFrame(iter_edge2ccNo(GG))
CCS.columns = ['sequence','modifier','accession','homo']
CCS = CCS.set_index(['sequence','modifier','accession'])
simpler_report['homo'] = CCS.homo



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

