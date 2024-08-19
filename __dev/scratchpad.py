import re
from math import inf

import duckdb
import furious_fastas as ff
import matplotlib.pyplot as plt
import networkx as nx
import numba
import numpy as np
import numpy.typing as npt
import pandas as pd
from networkx.algorithms import bipartite

fastas = ff.fastas(
    "data/Human_2024_02_16_UniProt_Taxon9606_Reviewed_20434entries_contaminant_tenzer.fasta"
)
fastas[0].sequence
duckcon = duckdb.connect()


# ms2rescore adapter
ms2rescore_input = pd.read_csv(
    "data/results.sage.ms2rescore.mokapot.peptides.txt", sep="\t"
)
ms2rescore_input = ms2rescore_input[
    [
        "peptide",
        "expmass",
        "retention_time",
        "charge",
        "mokapot q-value",
        "protein_list",
    ]
]


class SimpleReplace:
    def __init__(self, pattern: str = r"\[.*?\]"):
        self._pattern = re.compile(pattern)

    def apply(self, string: str, withwhat: str = ""):
        return re.sub(ptm_pattern, withwhat, string)


mods_bracket_anihilator = SimpleReplace()

ms2rescore_input["raw_sequence"] = ms2rescore_input.peptide.map(
    mods_bracket_anihilator.apply
)
peptide_report = ms2rescore_input

raw_sequences = set(ms2rescore_input["raw_sequence"])


@numba.njit(parallel=True)
def get_adjacency_matrix(strings: list[str], substrings: list[str]) -> npt.NDArray:
    res = np.empty(shape=(len(strings), len(substrings)), dtype=np.bool_)
    for string_id in numba.prange(len(strings)):
        for substring_id, substring in enumerate(substrings):
            res[string_id, substring_id] = substring in strings[string_id]
    return res


def get_edges(fastas: list[str], peptides: list[str]) -> list[tuple[int, int]]:
    adjacency_matrix = get_adjacency_matrix(fastas, peptides)
    res = []
    for prot_id, pep_id in zip(*np.nonzero(adjacency_matrix)):
        res.append((int(prot_id), int(pep_id)))
    res.sort()
    del adjacency_matrix
    return res


edges = get_edges([f.sequence for f in fastas], list(raw_sequences))


def invert_peptide_indices(edges):
    """
    Stupid networkx won't distinguish two sets of nodes.
    So we take peptide_id -> -1-peptide_id
    """
    return [(prot_id, -1 - pep_id) for prot_id, pep_id in edges]


edges = invert_peptide_indices(edges)

prots = set()
peps = set()
for prot, pep in edges:
    prots.add(prot)
    peps.add(pep)

pep2prot = nx.Graph()
pep2prot.add_nodes_from(list(prots), bipartite=0)
pep2prot.add_nodes_from(list(peps), bipartite=1)
pep2prot.add_edges_from(edges)

len(pep2prot.edges)
len(pep2prot)


for cc in nx.connected_components(pep2prot):
    if len(cc) > 100:
        break


def graph_is_lexicographically_sorted(edges):
    a_prev = -inf
    b_prev = -inf
    for a, b in edges:
        if a < a_prev or (a == a_prev and b < b_prev):
            return False
        a_prev = a
        b_prev = b
    return True


for cc in nx.connected_components(pep2prot):
    if not graph_is_lexicographically_sorted(pep2prot.subgraph(cc).edges):
        break

pep2prot.subgraph(cc).edges

all(
    graph_is_lexicographically_sorted(pep2prot.subgraph(cc).edges)
    for cc in nx.connected_components(pep2prot)
)  # weird...


all(map(graph_is_lexicographically_sorted, nx.connected_components(pep2prot)))


# Create a subgraph SG based on a (possibly multigraph) G


def induce(G, cc):
    SG = G.__class__()
    SG.add_nodes_from((n, G.nodes[n]) for n in cc)
    if SG.is_multigraph():
        SG.add_edges_from(
            (n, nbr, key, d)
            for n, nbrs in G.adj.items()
            if n in cc
            for nbr, keydict in nbrs.items()
            if nbr in cc
            for key, d in keydict.items()
        )
    else:
        SG.add_edges_from(
            (n, nbr, d)
            for n, nbrs in G.adj.items()
            if n in cc
            for nbr, d in nbrs.items()
            if nbr in cc
        )
    SG.graph.update(G.graph)
    return SG


induce(pep2prot, cc).edges
CC.adj
CC = pep2prot.subgraph(cc)


# CC = nx.subgraph(pep2prot, cc)
# trimming nodes where proteins are not good

# len(pep2prot.edges)
# nx.draw(CC, with_labels=True)
# plt.show()

# CC.edges


## VERY SLOW
# substrings = re.compile(
#     "|".join(set(ms2rescore_input.peptide.map(mods_bracket_anihilator.apply)))
# )
# big_fasta = " ".join(fasta.sequence for fasta in fastas)
# matches = re.finditer(substrings, big_fasta)

# # Iterate through matches and print their start and end positions
# for match in matches:
#     start = match.start()
#     end = match.end()
#     substring = match.group()


## SLOW AND NOT TRIVIAL TO PARALLELIZE IN NUMBA
# @numba.njit
# def test_quad(strings, substrings):
#     for string_id, string in enumerate(strings):
#         for substring_id, substring in enumerate(substrings):
#             if substring in string:
#                 res.append((string_id, substring_id))
#     return res


# test_quad([f.sequence for f in fastas], list(raw_sequences))
