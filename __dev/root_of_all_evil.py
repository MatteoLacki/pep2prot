%load_ext autoreload
%autoreload 2
import multiprocessing as mp
import typing
from collections import Counter, defaultdict

import networkx as nx
import numba
import numpy as np
import numpy.typing as npt
import pandas as pd
from numba import types
from numba.typed import Dict
from numba_progress import ProgressBar
from tqdm import tqdm

import ahocorasick
import furious_fastas as ff
from pep2prot.graph_ops import (asserting_unique, get_adjacency_matrix,
                                get_minimal_protein_group_coverage,
                                get_protein_group_cover_greadily,
                                get_protein_groups)
from pep2prot.readers import read_DiaNN

pd.set_option('display.max_columns', None)
fastas_path = "data/betterwheat/20211231_UniProtKB2021_04_4565TriticumAestivum_143648entries_172contaminants.fasta"

fastas = [(f.header.split(" ", 1)[0][1:], f.sequence) for f in ff.fastas(fastas_path)]
diann_report_extract = pd.read_parquet("data/diann_report_extract.parquet")

# covering_protein_groups, adjacency_matrix = get_minimal_protein_group_coverage(
#     peptides=list(set(diann_report_extract.peptide)),
#     fastas=fastas,
#     min_number_of_peptides=3,
# )


@numba.njit(boundscheck=True)
def make_fasta_mask(prot_lens: list[int]) -> npt.NDArray:
    fastas_mask = np.full(fill_value=-1, shape=sum(prot_lens) + len(prot_lens) - 1)
    i = 0
    for prot_id, prot_len in enumerate(prot_lens):
        fastas_mask[i:i+prot_len] = prot_id
        i += prot_len + 1
    return fastas_mask

def is_nondecreasing(xx: typing.Iterable) -> bool:
    xx = iter(xx)
    x_prev = next(xx)
    for x in xx:
        if x < x_prev:
            return False
        x_prev = x
    return True

def iter_prot_pep_edges(big_fasta: str, peptides: list[str], prot_ids_in_big_fasta: npt.NDArray,  check_ok: bool=True,) -> typing.Iterator[tuple[int,int]]:
    """Iterate over the edge of the protein-peptide graph."""
    automaton = ahocorasick.Automaton()
    for pep_id,pep in enumerate(peptides):
        automaton.add_word(pep, pep_id)
    automaton.make_automaton()
    for end, pep_id in automaton.iter(big_fasta):
        pep = peptides[pep_id]
        start = end - len(pep) + 1
        if check_ok:
            assert big_fasta[start:end+1] == pep, "ahocorasick.Automaton did not find the peptide."
        prot_id = prot_ids_in_big_fasta[start]
        yield int(prot_id), pep_id

def make_protein_peptide_graph(peptides: list[str], protein_sequences: list[str]) -> dict[int,int]:
    prot_ids_in_big_fasta = make_fasta_mask(prot_lens = list(map(len, protein_sequences)))
    big_fasta = " ".join(protein_sequences)
    graph = defaultdict(list)
    for prot_id, pep_id in iter_prot_pep_edges(big_fasta, peptides, prot_ids_in_big_fasta):
        neg_pep_id = -pep_id-1
        graph[prot_id].append(neg_pep_id)
        graph[neg_pep_id].append(prot_id)
    return dict(graph)

def fill_peptide_counts(proteins_peptide_cnt: npt.NDArray, graph: dict[int,int]) -> None:
    for prot_id, pep_ids in graph.items():
        if prot_id >= 0:
            proteins_peptide_cnt[prot_id] += len(pep_ids)

def fill_protein_groups(prot_to_protgroup: npt.NDArray, starting_group: int = -1) -> None:
    pep_ids_to_protgroup = {}
    group = starting_group
    for prot_id, pep_ids in graph.items():
        if prot_id >= 0:# proteins nonzero, peptides negative
            pep_ids = tuple(pep_ids)
            if not pep_ids in pep_ids_to_protgroup:
                group += 1
                pep_ids_to_protgroup[tuple(pep_ids)] = group
                prot_to_protgroup[prot_id] = group
            else:
                prot_to_protgroup[prot_id] = pep_ids_to_protgroup[pep_ids]

def get_protein_groups(proteins_with_enough_peptides: pd.DataFrame) -> pd.DataFrame:
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


def iter_pgrepr_pep_edges(graph: dict[int,int], protein_groups_representitives: typing.Iterable[int]) -> typing.Iterator[tuple[int,int]]:
    protein_groups_representitives = set(protein_groups_representitives)
    for prot_id, pep_ids in graph.items():
        if prot_id >= 0:
            for pep_id in pep_ids:
                yield prot_id, pep_id

def get_protein_group_covers(
    edges: list[tuple[int, int]],
    cpu_cnt: int = mp.cpu_count(),
    _progressbar_msg: str = "Covering Peptides with Protein Groups",
) -> list[list[int]]:
    pep_prot_graph = nx.Graph(edges)
    pep_prot_subgraphs = [
        pep_prot_graph.subgraph(cc).copy()
        for cc in nx.connected_components(pep_prot_graph)
    ]
    if _progressbar_msg != "":
        pep_prot_subgraphs = tqdm(pep_prot_subgraphs, desc=_progressbar_msg)
    with mp.Pool(cpu_cnt) as pool:
        covers = list(pool.map(get_protein_group_cover_greadily, pep_prot_subgraphs))

    return covers


# %%time
peptides = list(set(diann_report_extract.peptide))
min_number_of_peptides: int = 3
cpu_cnt: int = mp.cpu_count()

proteins = pd.DataFrame({"header": [header for header, sequence in fastas]})
assert len(set(proteins.header)) == len(proteins), "Protein headers are not unique."
proteins.index.name = "prot_id"

graph = make_protein_peptide_graph(peptides, protein_sequences=[seq for _, seq in fastas])
proteins["peptide_cnt"] = 0
fill_peptide_counts(proteins.peptide_cnt.to_numpy(), graph)

proteins["group"] = -1
fill_protein_groups(proteins.group.to_numpy())

proteins_with_enough_peptides = proteins[
    proteins.peptide_cnt >= min_number_of_peptides
]
protein_groups = get_protein_groups(proteins_with_enough_peptides)

# this could also factorize through connected components, so as to do multiple permutations and different covers.
protein_representatives_covers = get_protein_group_covers(edges=list(iter_pgrepr_pep_edges(graph, protein_groups.repr_prot_id)), cpu_cnt=cpu_cnt,
        _progressbar_msg="Covering Peptides with Protein Groups",)

protein_representatives_cover = {
    prot_id for cover in protein_representatives_covers for prot_id in cover
}
protein_groups = protein_groups[# retaining only PGs in the cover
    protein_groups.repr_prot_id.isin(protein_representatives_cover)
].copy()

protein_groups["pep_ids"] = [
    [-inv_pep_id-1 for inv_pep_id in graph[prot_id]]
    for prot_id in protein_groups.repr_prot_id
]

assert all(
    pep_cnt == len(pep_ids)
    for pep_cnt, pep_ids in zip(protein_groups.peptide_cnt, protein_groups.pep_ids)
), "Some peptide groups had different number of peptides than anticipated."

initial_covered_pep_ids = set(pep_id for prot_id, pep_ids in graph.items() if prot_id < 0 for pep_id in pep_ids)
len(initial_covered_pep_ids)

final_covered_pep_ids = set(pep_id for pep_ids in protein_groups.pep_ids for pep_id in pep_ids)
assert initial_covered_pep_ids == final_covered_pep_ids, "Peptides covered by proteins do not contain the same peptides as those induced by protein cover."
