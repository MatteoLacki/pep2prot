import math
import multiprocessing as mp
import typing
from collections import Counter, defaultdict
from math import inf

import networkx as nx
import numba
import numpy as np
import numpy.typing as npt
import pandas as pd
from numba_progress import ProgressBar
from tqdm import tqdm

import ahocorasick
import furious_fastas as ff
from pep2prot.misc import asserting_unique, is_nondecreasing
from pep2prot.numpy_ops import check_arrays_equal, check_vecs_the_same


@numba.njit(parallel=True)
def get_adjacency_matrix(
    strings: list[str],
    substrings: list[str],
    progres_proxy: ProgressBar,
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
            progres_proxy.update(1)
    return res


@numba.njit
def faster_nonzero(adjacency_matrix):
    return np.nonzero(adjacency_matrix)


@numba.njit
def sparsify_adjacency_matrix(
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


def get_protein_group_cover_greedily(bipartite_graph: nx.Graph) -> list[int]:
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


def get_protein_group_covers(
    edges: list[tuple[int, int]],
    cpu_cnt: int = mp.cpu_count(),
    _progressbar_msg: str = "Covering Peptides with Protein Groups",
) -> list[list[int]]:
    # only need to change peps to -1-pep cause inverted peptides can be discarded after properly named protein groups forming a cover are selected
    pep_prot_graph = nx.Graph(invert_peptide_indices(edges))
    pep_prot_subgraphs = [
        pep_prot_graph.subgraph(cc).copy()
        for cc in nx.connected_components(pep_prot_graph)
    ]
    if _progressbar_msg != "":
        pep_prot_subgraphs = tqdm(pep_prot_subgraphs, desc=_progressbar_msg)
    with mp.Pool(cpu_cnt) as pool:
        covers = list(pool.map(get_protein_group_cover_greedily, pep_prot_subgraphs))

    return covers


def get_minimal_protein_group_coverage_slow(
    peptides: list[str],
    fastas: list[tuple[str, str]],
    min_number_of_peptides: int = 3,
    cpu_cnt: int = mp.cpu_count(),
) -> tuple[pd.DataFrame, npt.NDArray]:
    """
    Get the approximate minimal protein group cover of submitted peptides slowly.

    Explicitly constructs the adjacency matrix.

    Arguments:
        peptide_report (list[str]): A list of found peptide sequences.
        fastas (list[tuple[str,str]]): A list of tuples with header and protein sequence each.
        min_number_of_peptides (int): Minimal number of peptides per protein to be even worthy of being considered a cover.
        cpu_cnt (int): Number of workers in a multiprocessing Pool to be used.

    Returns:
        tuple[pd.DataFrame, npt.NDArray]: A table containing protein groups in the greedily approximated minimal protein group cover and the adjacency matrix.
    """
    protein_sequences = [sequence for header, sequence in fastas]
    proteins = pd.DataFrame({"header": [header for header, sequence in fastas]})
    assert len(set(proteins.header)) == len(proteins), "Protein headers are not unique."
    proteins.index.name = "prot_id"
    with ProgressBar(
        total=len(proteins) * len(peptides), desc="Getting Adjacency Matrix"
    ) as progres_proxy:
        adjacency_matrix = get_adjacency_matrix(
            protein_sequences, peptides, progres_proxy
        )

    proteins["peptide_cnt"] = adjacency_matrix.sum(axis=1)
    # optimization: first filter by peptide counts first.
    proteins["group"] = get_protein_groups(np.packbits(adjacency_matrix, axis=-1))
    assert check_protein_groups_are_OK(
        proteins.group.to_numpy(),
        adjacency_matrix,
    ), "Some proteins in the same protein group do not share the same peptides (pattern of 0s and 1s differs)."
    proteins_with_enough_peptides = proteins[
        proteins.peptide_cnt >= min_number_of_peptides
    ]
    protein_groups = (
        proteins_with_enough_peptides.reset_index()
        .groupby("group")
        .aggregate(
            {
                "prot_id": [
                    list,
                    "min",
                ],
                "peptide_cnt": asserting_unique,
                "header": list,
                "group": len,
            }
        )
    )
    protein_groups.columns = [
        "prot_idxs",
        "repr_prot_id",
        "peptide_cnt",
        "accessions",
        "prot_in_group_cnt",
    ]

    PG_adjmat = adjacency_matrix[protein_groups.repr_prot_id.to_numpy()]
    assert np.all(
        protein_groups.peptide_cnt == PG_adjmat.sum(axis=1)
    ), "Quality check failed."

    prot_pep_edges = sparsify_adjacency_matrix(
        PG_adjmat, protein_groups.repr_prot_id.to_numpy()
    )
    protein_representatives_covers = get_protein_group_covers(
        prot_pep_edges,
        cpu_cnt=cpu_cnt,
        _progressbar_msg="Covering Peptides with Protein Groups",
    )

    PG_covering_connected_component_dist = Counter(
        map(len, protein_representatives_covers)
    )
    protein_representatives_cover = {
        prot_id for cover in protein_representatives_covers for prot_id in cover
    }
    protein_groups = protein_groups[
        protein_groups.repr_prot_id.isin(protein_representatives_cover)
    ].copy()

    protein_groups["pep_ids"] = protein_groups.repr_prot_id.map(
        pd.DataFrame(sorted(prot_pep_edges), columns=["prot_id", "pep_id"])
        .groupby("prot_id")
        .aggregate(dict(pep_id=list))
        .pep_id
    )

    assert all(
        pep_cnt == len(pep_ids)
        for pep_cnt, pep_ids in zip(protein_groups.peptide_cnt, protein_groups.pep_ids)
    ), "Some peptide groups had different number of peptides than anticipated."

    return protein_groups, adjacency_matrix


@numba.njit(boundscheck=True)
def make_fasta_mask(prot_lens: list[int]) -> npt.NDArray:
    """
    Make an array mapping each letter in a big concatenated fasta string to it fasta number in the order of appearance in the file.
    """
    fastas_mask = np.full(fill_value=-1, shape=sum(prot_lens) + len(prot_lens) - 1)
    i = 0
    for prot_id, prot_len in enumerate(prot_lens):
        fastas_mask[i : i + prot_len] = prot_id
        i += prot_len + 1
    return fastas_mask


def iter_prot_pep_edges(
    big_fasta: str,
    peptides: list[str],
    prot_ids_in_big_fasta: npt.NDArray,
    check_ok: bool = True,
) -> typing.Iterator[tuple[int, int]]:
    """
    Iterate over the edges of the protein-peptide graph.

    Arguments:
        big_fasta (str): concatenated protein sequences separated by space.
        peptides (list[str]): Observed peptide strings.
        prot_ids_in_big_fasta (npt.NDArray): an array mapping each letter in a big concatenated fasta string to it fasta number in the order of appearance in the file.
        check_ok (bool): Check the ahocorasick Automaton produces valid subsequences.

    Yield:
        tuple[int,int]: Protein and peptide id.
    """
    automaton = ahocorasick.Automaton()
    for pep_id, pep in enumerate(peptides):
        automaton.add_word(pep, pep_id)
    automaton.make_automaton()
    for end, pep_id in automaton.iter(big_fasta):
        pep = peptides[pep_id]
        start = end - len(pep) + 1
        if check_ok:
            assert (
                big_fasta[start : end + 1] == pep
            ), "ahocorasick.Automaton did not find the peptide."
        prot_id = prot_ids_in_big_fasta[start]
        yield int(prot_id), pep_id


def make_protein_peptide_graph(
    peptide_sequences: list[str],
    protein_sequences: list[str],
) -> dict[int, list[int]]:
    """Make a protein-peptide bipartite graph.

    Here one set of nodes consists of proteins, another of peptides, and edge corresponds to peptide sequence being a subsequence of a protein sequence.

    Arguments:
        peptide_sequences (list[str]): List of peptide sequences.
        protein_sequences (list[str]): List of protein sequences.

    Returns:
        dict: The mapping of nodes to other nodes in a graph. Strictly negative numbers correspond to peptide ids. To get a proper peptide id, follow 'pep_id = -inv_pep_id - 1'.
    """
    prot_ids_in_big_fasta = make_fasta_mask(prot_lens=list(map(len, protein_sequences)))
    big_fasta = " ".join(protein_sequences)
    graph = defaultdict(list)
    for prot_id, pep_id in iter_prot_pep_edges(
        big_fasta, peptide_sequences, prot_ids_in_big_fasta
    ):
        neg_pep_id = -pep_id - 1
        graph[prot_id].append(neg_pep_id)
        graph[neg_pep_id].append(prot_id)
    return dict(graph)


def fill_peptide_counts(
    proteins_peptide_cnt: npt.NDArray,
    graph: dict[int, list[int]],
) -> None:
    """Fill `proteins_peptide_cnt` with counts of peptides corresponding to protein ids in graph."""
    for prot_id, pep_ids in graph.items():
        if prot_id >= 0:
            proteins_peptide_cnt[prot_id] += len(pep_ids)


def fill_protein_groups(
    graph: dict[int, list[int]],
    prot_to_protgroup: npt.NDArray,
    starting_group: int = -1,
) -> None:
    pep_ids_to_protgroup: dict[tuple[int, ...], int] = {}
    group = starting_group
    for prot_id, pep_ids in graph.items():
        if prot_id >= 0:  # proteins nonzero, peptides negative
            pep_ids: tuple[int, ...] = tuple(pep_ids)
            if not pep_ids in pep_ids_to_protgroup:
                group += 1
                pep_ids_to_protgroup[tuple(pep_ids)] = group
                prot_to_protgroup[prot_id] = group
            else:
                prot_to_protgroup[prot_id] = pep_ids_to_protgroup[pep_ids]


def get_protein_groups_df(proteins_with_enough_peptides: pd.DataFrame) -> pd.DataFrame:
    """Aggregate protein info to protein groups summary."""
    protein_groups = (
        proteins_with_enough_peptides.reset_index()
        .groupby("group")
        .aggregate(
            {
                "prot_id": [
                    list,
                    "min",
                ],
                "peptide_cnt": asserting_unique,
                "header": list,
                "group": len,
            }
        )
    )
    protein_groups.columns = [
        "prot_idxs",
        "repr_prot_id",
        "peptide_cnt",
        "accessions",
        "prot_in_group_cnt",
    ]
    return protein_groups


def iter_pgrepr_pep_edges(
    graph: dict[int, list[int]],
    protein_groups_representitives: typing.Iterable[int],
) -> typing.Iterator[tuple[int, int]]:
    """Iterate over edges in the protein-peptide graph that correspond only to protein group representatives.

    Typically, we choose the first available protein in a group as representative.
    """
    protein_groups_representitives = set(protein_groups_representitives)
    for prot_id, pep_ids in graph.items():
        if prot_id >= 0 and prot_id in protein_groups_representitives:
            for pep_id in pep_ids:
                yield prot_id, pep_id


def get_protein_group_covers_2(
    edges: list[tuple[int, int]],
    cpu_cnt: int = mp.cpu_count(),
    _progressbar_msg: str = "Covering Peptides with Protein Groups",
) -> list[list[int]]:
    """
    Establish the minimal protein group cover of peptides covered by all protein groups.
    """
    pep_prot_graph = nx.Graph(edges)
    pep_prot_subgraphs = [
        pep_prot_graph.subgraph(cc).copy()
        for cc in nx.connected_components(pep_prot_graph)
    ]
    if _progressbar_msg != "":
        pep_prot_subgraphs = tqdm(pep_prot_subgraphs, desc=_progressbar_msg)
    with mp.Pool(cpu_cnt) as pool:
        covers = list(pool.map(get_protein_group_cover_greedily, pep_prot_subgraphs))

    return covers


def get_minimal_protein_group_coverage_fast(
    peptides: list[str],
    fastas: list[tuple[str, str]],
    min_number_of_peptides: int = 3,
    cpu_cnt: int = mp.cpu_count(),
    _debug: bool = False,
) -> (
    pd.DataFrame
    | tuple[pd.DataFrame, pd.DataFrame, dict[int, list[int]], list[list[int]]]
):
    """
    Get the approximate minimal protein group cover of submitted peptides.

    Using the Aho-Corasick algorithm underneath for peptide-to-protein graph build up.

    Arguments:
        peptide_report (list[str]): A list of found peptide sequences.
        fastas (list[tuple[str,str]]): A list of tuples with header and protein sequence each.
        min_number_of_peptides (int): Minimal number of peptides per protein to be even worthy of being considered a cover.
        cpu_cnt (int): Number of workers in a multiprocessing Pool to be used.
        _debug (bool): Return intermediate data structures.

    Returns:
        pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame, dict[int,list[int]], list[list[int]]] : A table containing protein groups in the greedily approximated minimal protein group cover. Optionally, that and a data frame with protein calculations results, the peptide-protein graph, and the list of covers per connected component.
    """

    proteins = pd.DataFrame({"header": [header for header, sequence in fastas]})
    proteins.index.name = "prot_id"

    assert len(set(proteins.header)) == len(proteins), "Protein headers are not unique."

    graph = make_protein_peptide_graph(
        peptides, protein_sequences=[seq for _, seq in fastas]
    )

    proteins["peptide_cnt"] = 0
    fill_peptide_counts(proteins.peptide_cnt.to_numpy(), graph)

    proteins[
        "group"
    ] = (
        starting_group
    ) = -1  # sentinel: -1 means a protein was not matched with any peptide
    fill_protein_groups(graph, proteins.group.to_numpy(), starting_group=starting_group)

    proteins_with_enough_peptides = proteins[
        proteins.peptide_cnt >= min_number_of_peptides
    ]

    protein_groups = get_protein_groups_df(proteins_with_enough_peptides)

    # this could also factorize through connected components, so as to do multiple permutations and different covers.
    representative_edges = list(
        iter_pgrepr_pep_edges(graph, protein_groups.repr_prot_id)
    )

    protein_representatives_covers = get_protein_group_covers_2(
        edges=representative_edges,
        cpu_cnt=cpu_cnt,
        _progressbar_msg="Covering Peptides with Protein Groups",
    )

    protein_representatives_cover = {
        prot_id for cover in protein_representatives_covers for prot_id in cover
    }

    protein_groups = protein_groups[  # retaining only PGs in the cover
        protein_groups.repr_prot_id.isin(protein_representatives_cover)
    ].copy()

    protein_groups["pep_ids"] = [
        [-inv_pep_id - 1 for inv_pep_id in graph[prot_id]]
        for prot_id in protein_groups.repr_prot_id
    ]

    assert all(
        pep_cnt == len(pep_ids)
        for pep_cnt, pep_ids in zip(protein_groups.peptide_cnt, protein_groups.pep_ids)
    ), "Some peptide groups had different number of peptides than anticipated."

    pep_ids_before_cover = set(
        -inv_pep_id - 1 for prot_id, inv_pep_id in representative_edges
    )
    final_covered_pep_ids = set(
        pep_id for pep_ids in protein_groups.pep_ids for pep_id in pep_ids
    )
    assert (
        pep_ids_before_cover == final_covered_pep_ids
    ), "Protein covering algorithm lost some peptides."

    if _debug:
        return protein_groups, proteins, graph, protein_representatives_covers

    return protein_groups
