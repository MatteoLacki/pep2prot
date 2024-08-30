%load_ext autoreload
%autoreload 2
import multiprocessing as mp
import typing
from collections import Counter

import networkx as nx
import numba
import numpy as np
import numpy.typing as npt
import pandas as pd
from tqdm import tqdm

import furious_fastas as ff
from pep2prot.graph_ops import (get_minimal_protein_group_coverage_fast,
                                get_protein_group_cover_greedily,
                                get_protein_groups_df, iter_prot_pep_edges,
                                make_fasta_mask, make_protein_peptide_graph,
                                max_degree_node)

pd.set_option('display.max_columns', None)
fastas_path = "data/betterwheat/20211231_UniProtKB2021_04_4565TriticumAestivum_143648entries_172contaminants.fasta"

fastas = [(f.header.split(" ", 1)[0][1:], f.sequence) for f in ff.fastas(fastas_path)]
diann_report_extract = pd.read_parquet("data/diann_report_extract.parquet")

# %%timeit
# protein_groups, proteins, graph, protein_representatives_covers = get_minimal_protein_group_coverage_fast(
#     peptides = list(set(diann_report_extract.peptide)),
#     fastas = fastas,
#     _debug = True,
# )


def make_protein_peptide_graph(
    peptide_sequences: list[str],
    protein_sequences: list[str],
) -> nx.Graph:
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
    graph = nx.Graph()
    for prot_id, pep_id in iter_prot_pep_edges(
        big_fasta, peptide_sequences, prot_ids_in_big_fasta
    ):
        neg_pep_id = -pep_id - 1
        graph.add_edge(prot_id, neg_pep_id)
    return graph

def fill_peptide_counts(
    proteins_peptide_cnt: npt.NDArray[np.intp],
    graph: nx.Graph,
) -> None:
    """Fill proteins_peptide_cnt` with counts of peptides corresponding to protein ids in graph."""
    for edge in graph.edges:
        prot_id, pep_id = edge if edge[0] >= 0 else edge[::-1]
        proteins_peptide_cnt[prot_id] += 1

def fill_protein_groups(
    graph: nx.Graph,
    prot_to_protgroup: npt.NDArray,
    starting_group: int = -1,
) -> None:
    """For each protein, write down its protein group number.

    The protein group numbers are defined simply as consecutive different tuples of 
    """
    pep_ids_to_protgroup: dict[tuple[int, ...], int] = {}
    group = starting_group
    for edge in graph.edges():
        prot_id, pep_id = edge if edge[0] >= 0 else edge[::-1]
        pep_ids: tuple[int, ...] = tuple(sorted(graph[prot_id])) # order of peptides irrelevant
        if not pep_ids in pep_ids_to_protgroup:
            group += 1
            pep_ids_to_protgroup[tuple(pep_ids)] = group
            prot_to_protgroup[prot_id] = group
        else:
            prot_to_protgroup[prot_id] = pep_ids_to_protgroup[pep_ids]


def iter_pgrepr_pep_edges(
    graph: nx.Graph,
    protein_groups_representitives: pd.Series | npt.NDArray[np.intp],
) -> typing.Iterator[tuple[int, int]]:
    """Iterate over edges in the protein-peptide graph that correspond only to protein group representatives.

    Typically, we choose the first available protein in a group as representative.
    """
    protein_groups_representitives = set(protein_groups_representitives)
    for edge in graph.edges():
        prot_id, pep_id = edge if edge[0] >= 0 else edge[::-1]
        if prot_id in protein_groups_representitives:
            for pep_id in graph[prot_id]:
                yield prot_id, pep_id


def find_connected_components(graph: nx.Graph, verbose: bool= True):
    ccs = nx.connected_components(graph)
    if verbose:
        ccs = tqdm(ccs, desc="Finding connected components.")
    return [ graph.subgraph(cc).copy() for cc in ccs ]

@numba.njit
def harmonic_number(n: int) -> float:
    assert n > 0, "Harmonic numbers are only defined for positive numbers."
    res = 0
    for k in range(1,n+1):
        res += 1/k
    return res

def get_protein_group_cover_greedily(bipartite_graph: nx.Graph) -> tuple[list[int], float]:
    cover = []
    peptides_to_cover = sum(n < 0 for n in bipartite_graph.nodes)
    max_cardinality_set_size = 0
    while peptides_to_cover > 0:
        prot_id, degree = max_degree_node(bipartite_graph)
        max_cardinality_set_size = max(max_cardinality_set_size, degree)
        cover.append(prot_id)
        peptides_to_cover -= degree
        for pep_id in list(bipartite_graph[prot_id]):
            bipartite_graph.remove_node(pep_id)
        bipartite_graph.remove_node(prot_id)

    assert len(bipartite_graph.edges) == 0, "Graph still has some edges which is fishy."
    return cover, max_cardinality_set_size

def get_covers(pep_protgroup_connected_components: list[nx.Graph]) -> pd.DataFrame:
    _ = harmonic_number(10)# just heating up numba
    with mp.Pool(cpu_cnt) as pool:
        covers = pd.DataFrame(
            pool.map(
                get_protein_group_cover_greedily,
                pep_protgroup_connected_components,
            ),
            columns=["cover", "max_cardinality_set_size"]
        )
    covers["cover_size"] = covers.cover.map(len)
    covers["harmonic_number"] = covers.cover_size.map(harmonic_number)
    covers["lowest_minimal_cover_size"] = np.ceil( covers.cover_size / covers.harmonic_number).astype(int)
    return covers


def get_minimal_protein_group_coverage(
    peptides: list[str],
    fastas: list[tuple[str, str]],
    min_number_of_peptides: int = 3,
    cpu_cnt: int = mp.cpu_count(),
    _debug: bool = False,
) -> (
    pd.DataFrame
    | tuple[pd.DataFrame, pd.DataFrame, nx.Graph, nx.Graph, pd.DataFrame, pd.DataFrame] 
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
        pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame, nx.Graph, nx.Graph, pd.DataFrame, pd.DataFrame]: A table containing protein groups in the greedily approximated minimal protein group cover. Optionally, that and a data frame with protein calculations results, the peptide-protein graph, the induce protein representative-peptide graph, a data frame with cover statistics, and a peptide view data frame.
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
    ) = -1
    # sentinel: -1 means a protein was not matched with any peptide
    fill_protein_groups(graph, proteins.group.to_numpy(), starting_group=starting_group)

    proteins_with_enough_peptides = proteins[
        proteins.peptide_cnt >= min_number_of_peptides
    ]

    protein_groups = get_protein_groups_df(proteins_with_enough_peptides)

    pep_protgroup_graph = nx.Graph(
        iter_pgrepr_pep_edges(graph, protein_groups.repr_prot_id)
    )
    pep_protgroup_connected_components = find_connected_components(pep_protgroup_graph)

    covers = get_covers(pep_protgroup_connected_components)
    protein_representatives_cover = {
        prot_id for cover in covers.cover for prot_id in cover
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

    pep_ids_before_cover = set(-n-1 for n in pep_protgroup_graph.nodes if n < 0)
    final_covered_pep_ids = set(
        pep_id for pep_ids in protein_groups.pep_ids for pep_id in pep_ids
    )
    assert (
        pep_ids_before_cover == final_covered_pep_ids
    ), "Protein covering algorithm lost some peptides."

    if _debug:
        def iter_edges_between_representatives_and_pep_ids(pep_protgroup_graph, protein_representatives_cover):
            for edge in pep_protgroup_graph.edges:
                prot_id, pep_id = edge if edge[0] >= 0 else edge[::-1]
                if prot_id in protein_representatives_cover:
                    yield prot_id, -pep_id-1

        final_edges = pd.DataFrame(iter_edges_between_representatives_and_pep_ids(pep_protgroup_graph, protein_representatives_cover), columns=["prot_id","pep_id"])

        repr_prot_id_to_protein_group = {
            repr_prot_id: prot_id for repr_prot_id, prot_id in zip(protein_groups.repr_prot_id, protein_groups.index)
        }
        final_edges["protein_group"] = final_edges.prot_id.map(repr_prot_id_to_protein_group)

        peptide_view = final_edges.groupby("pep_id").agg({"protein_group":list})
        peptide_view["protein_group_size"] = peptide_view.protein_group.map(len)

        return protein_groups, proteins, graph, pep_protgroup_graph, covers, peptide_view
    return protein_groups

peptides = list(set(diann_report_extract.peptide))
