import math
import multiprocessing as mp
import re
from collections import Counter
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
        return re.sub(self._pattern, withwhat, string)


x = len(substrings) // 64


@numba.njit(parallel=True)
def get_adjacency_matrix(
    strings: list[str], substrings: list[str], _padding: int = 64
) -> npt.NDArray:
    """
    Get the adjacency matrix for the problem.
    It is padded to the right with 0s.
    """
    assert _padding in (32, 64), "Our code works on 32 and 64 bit architectures."
    div, mod = np.divmod(len(substrings), _padding)
    res = np.zeros(
        shape=(
            len(strings),
            (div + (mod > 0)) * _padding,
        ),
        dtype=np.bool_,
    )
    for string_id in numba.prange(len(strings)):
        for substring_id, substring in enumerate(substrings):
            res[string_id, substring_id] = substring in strings[string_id]
    return res


@numba.njit
def faster_nonzero(adjacency_matrix):
    return np.nonzero(adjacency_matrix)


@numba.njit
def sparsify_adjacency_martrix(
    adjacency_matrix: npt.NDArray,
    protein_groups_representatives: npt.NDArray,
) -> list[tuple[int, int]]:
    return [
        (protein_groups_representatives[i], j)
        for i, j in zip(*faster_nonzero(adjacency_matrix))
    ]


@numba.njit
def invert_peptide_indices(edges):
    """
    Stupid networkx won't distinguish two sets of nodes.
    So we take peptide_id -> -1-peptide_id
    """
    return [(prot_id, -1 - pep_id) for prot_id, pep_id in edges]


def edges_to_bipartite_graph(edges: list[tuple[int, int]]) -> nx.Graph:
    prots = set()
    peps = set()
    for prot, pep in edges:
        prots.add(prot)
        peps.add(pep)
    G = nx.Graph()
    G.add_nodes_from(list(prots), bipartite=0)
    G.add_nodes_from(list(peps), bipartite=1)
    G.add_edges_from(edges)
    return G


@numba.njit
def graph_nodes_to_proteins_and_peptides_indices_lists(
    graph_nodes,
) -> tuple[list[int], list[int]]:
    prot_rows = []
    pep_cols = []
    for node in graph_nodes:
        if node < 0:
            pep_cols.append(-node - 1)
        else:
            prot_rows.append(node)
    pep_cols.sort()
    prot_rows.sort()
    return prot_rows, pep_cols


# not necessary
@numba.njit
def get_subadjacency_matrix(
    adjacency_matrix: npt.NDArray, prot_rows: list[int], pep_cols: list[int]
) -> npt.NDArray:
    res = np.zeros(shape=(len(prot_rows), len(pep_cols)), dtype=np.bool_)
    for i, prot_id in enumerate(prot_rows):
        for j, pep_id in enumerate(pep_cols):
            res[i, j] = adjacency_matrix[prot_id, pep_id]
    return res


# it seems we cannot do it in numba this way.
def get_row_sorting_order(matrix):
    assert len(matrix.shape) == 2, "Wrong dimension: expecting 2D numpy array."
    return np.argsort(
        matrix.view([("", matrix.dtype)] * matrix.shape[1]), axis=0
    ).ravel()


@numba.njit
def check_arrays_equal(a, b):
    assert len(a) == len(b)
    for i in range(len(a)):
        if a[i] != b[i]:
            return False
    return True


@numba.njit
def group_rows(rows, row_order):
    """Collect rows that are the same in the lexicographical order.

    Returns:
        np.array: An assignment of row_order indices into groups starting from 0 onwards.
    """
    groups = np.empty(shape=row_order.shape, dtype=np.uint32)
    groups[0] = 0
    current_group = 0
    for i in range(1, len(row_order)):
        current_group += not check_arrays_equal(
            rows[row_order[i]],
            rows[row_order[i - 1]],
        )
        groups[i] = current_group
    return groups


def get_dtype_and_its_bits():
    for x in (32, 64):
        if str(x) in str(np.uintp):
            return np.uintp, x


@numba.njit
def simplify_groups(xx, yy):
    mapping = {}
    result = np.empty(dtype=np.uint32, shape=len(xx))
    i = 0
    for j, (x, y) in enumerate(zip(xx, yy)):
        if (x, y) not in mapping:
            mapping[(x, y)] = i
            i += 1
        result[j] = mapping[(x, y)]
    return result


@numba.njit
def get_groups(group_counters):
    partial_groups = group_counters[:, 0].copy().astype(np.uint32)
    for j in range(1, group_counters.shape[1]):
        partial_groups = simplify_groups(partial_groups, group_counters[:, j])
    return partial_groups


@numba.njit
def check_vecs_the_same(array):
    if len(array) <= 1:
        return True
    prev = array[0]
    for i in range(len(array)):
        if not np.all(prev == array[i]):
            return False
        prev = array[i]
    return True


@numba.njit
def check_groups_are_OK(protein_groups, adjacency_matrix):
    """For testing if the groups corresponding to an adjacency matrix are good."""
    for i in range(max(protein_groups)):
        W = adjacency_matrix[protein_groups == i]
        if not check_vecs_the_same(W):
            return False
    return True


def max_degree_node(bipartite_graph) -> tuple[int, int]:
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


def graph_is_lexicographically_sorted(edges):
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

mods_bracket_anihilator = SimpleReplace()
ms2rescore_input["raw_sequence"] = ms2rescore_input.peptide.map(
    mods_bracket_anihilator.apply
)
peptide_report = ms2rescore_input
raw_sequences = list(set(ms2rescore_input["raw_sequence"]))
adjacency_matrix = get_adjacency_matrix([f.sequence for f in fastas], raw_sequences)

# get the protein group indices:
peptides_per_protein_cnt = adjacency_matrix.sum(axis=1)
protein_groups = get_groups(np.packbits(adjacency_matrix, axis=-1))
assert check_groups_are_OK(protein_groups, adjacency_matrix)


# filtering graph: can do that before other calculations? no.
proteins_with_enough_peptides = np.nonzero(
    peptides_per_protein_cnt >= min_number_of_peptides
)[0]
protein_groups_representatives = np.sort(np.unique(proteins_with_enough_peptides))
reps_adjacency_matrix = adjacency_matrix[protein_groups_representatives]

assert np.all(
    peptides_per_protein_cnt[protein_groups_representatives]
    == reps_adjacency_matrix.sum(axis=1)
), "Quality check failed."

edges = sparsify_adjacency_martrix(
    reps_adjacency_matrix, protein_groups_representatives
)


get_protein_group_cover(edges)


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
