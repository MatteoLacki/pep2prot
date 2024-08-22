import math
import multiprocessing as mp
from collections import Counter
from math import inf

import furious_fastas as ff
import matplotlib.pyplot as plt
import networkx as nx
import numba
import numpy as np
import numpy.typing as npt
import pandas as pd
from networkx.algorithms import bipartite

from pep2prot.graph_ops import get_adjacency_matrix
from pep2prot.readers import SimpleReplace


def max_degree_node(bipartite_graph: nx.Graph) -> tuple[int, int]:
    if len(bipartite_graph) == 0:
        return 0
    degree = -1
    most_covering_node = -1
    for node in bipartite_graph.nodes:
        if node >= 0:
            node_deg = nx.degree(bipartite_graph, node)
            if node_deg >= degree:
                degree = node_deg
                most_covering_node = node
    assert most_covering_node != -1, "There were no proteins in the graph."
    return most_covering_node, degree


def get_protein_group_cover_greadily(bipartite_graph: nx.Graph) -> list[int]:
    cover = []
    peptides_to_cover = sum(n < 0 for n in bipartite_graph.nodes)
    while peptides_to_cover > 0:
        prot_id, degree = max_degree_node(bipartite_graph)
        cover.append(prot_id)
        peptides_to_cover -= degree
        for pep_id in list(bipartite_graph[prot_id]):
            bipartite_graph.remove_node(pep_id)
        bipartite_graph.remove_node(prot_id)

    assert len(bipartite_graph.edges) == 0, "Graph still has some edges which is fishy."
    return cover


def graph_is_lexicographically_sorted(edges: list[int, int]) -> bool:
    a_prev = -inf
    b_prev = -inf
    for a, b in edges:
        if a < a_prev or (a == a_prev and b < b_prev):
            return False
        a_prev = a
        b_prev = b
    return True


def get_protein_group_cover(edges: list[tuple[int, int]]) -> list[int]:
    # only need to change peps to -1-pep cause inverted peptides can be discarded after properly named protein groups forming a cover are selected
    pep_prot_graph = nx.Graph(invert_peptide_indices(edges))
    pep_prot_subgraphs = [
        pep_prot_graph.subgraph(cc).copy()
        for cc in nx.connected_components(pep_prot_graph)
    ]
    with mp.Pool(cpu_cnt) as pool:
        covers = list(pool.map(get_protein_group_cover_greadily, pep_prot_subgraphs))

    cover = [prot_id for cover in covers for prot_id in cover]
    cover.sort()

    return cover


# THE AMOK OF SIMPLICITY
min_number_of_peptides = 3
cpu_cnt = mp.cpu_count()


fastas = ff.fastas(
    "data/Human_2024_02_16_UniProt_Taxon9606_Reviewed_20434entries_contaminant_tenzer.fasta"
)


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


mods_bracket_anihilator = SimpleReplace()
ms2rescore_input["raw_sequence"] = ms2rescore_input.peptide.map(
    mods_bracket_anihilator.apply
)
peptide_report = ms2rescore_input
raw_sequences = list(set(ms2rescore_input["raw_sequence"]))
adjacency_matrix = get_adjacency_matrix([f.sequence for f in fastas], raw_sequences)

peptides_per_protein_cnt = adjacency_matrix.sum(axis=1)
protein_groups = get_protein_groups(np.packbits(adjacency_matrix, axis=-1))
assert check_groups_are_OK(protein_groups, adjacency_matrix)

proteins_with_enough_peptides = np.nonzero(
    peptides_per_protein_cnt >= min_number_of_peptides
)[0]
protein_groups_representatives = np.sort(np.unique(proteins_with_enough_peptides))
reps_adjacency_matrix = adjacency_matrix[protein_groups_representatives]

assert np.all(
    peptides_per_protein_cnt[protein_groups_representatives]
    == reps_adjacency_matrix.sum(axis=1)
), "Quality check failed."

prot_pep_edges = sparsify_adjacency_martrix(
    reps_adjacency_matrix, protein_groups_representatives
)
protein_group_cover = get_protein_group_cover(prot_pep_edges)


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


# def test_packbits():
#     test = np.array([
#         [0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,1,1, 0,0,0,0,1,0,0,1,],
#         [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,1, 1,0,0,0,0,0,1,1, 1,0,0,0,1,0,0,1,],
#     ], dtype=np.bool_)
#     res = np.array([[  1,   0,   3,   9],
#            [  0,   1, 131, 137]], dtype=np.uint8)
#     assert np.all(np.packbits(test, axis=-1) == res)
