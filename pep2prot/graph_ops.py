import math
import multiprocessing as mp
from collections import Counter
from math import inf

import furious_fastas as ff
import networkx as nx
import numba
import numpy as np
import numpy.typing as npt

from pep2prot.numpy_ops import check_arrays_equal, check_vecs_the_same


@numba.njit(parallel=True)
def get_adjacency_matrix(
    strings: list[str],
    substrings: list[str],
    _padding: int = 64,
) -> npt.NDArray:
    """
    Get the adjacency matrix of the bipartite protein-peptide graph.
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


@numba.njit
def get_subadjacency_matrix(
    adjacency_matrix: npt.NDArray, prot_rows: list[int], pep_cols: list[int]
) -> npt.NDArray:
    res = np.zeros(shape=(len(prot_rows), len(pep_cols)), dtype=np.bool_)
    for i, prot_id in enumerate(prot_rows):
        for j, pep_id in enumerate(pep_cols):
            res[i, j] = adjacency_matrix[prot_id, pep_id]
    return res


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


@numba.njit
def _simplify_groups(xx, yy):
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
def get_protein_groups(group_counters: npt.NDArray):
    partial_groups = group_counters[:, 0].copy().astype(np.uint32)
    for j in range(1, group_counters.shape[1]):
        partial_groups = _simplify_groups(partial_groups, group_counters[:, j])
    return partial_groups


@numba.njit
def check_protein_groups_are_OK(protein_groups, adjacency_matrix):
    """For testing if the groups corresponding to an adjacency matrix are good."""
    for i in range(max(protein_groups)):
        W = adjacency_matrix[protein_groups == i]
        if not check_vecs_the_same(W):
            return False
    return True


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


def get_protein_group_cover(
    edges: list[tuple[int, int]],
    cpu_cnt: int = mp.cpu_count(),
) -> list[int]:
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
