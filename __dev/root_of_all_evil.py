%load_ext autoreload
%autoreload 2
import pandas as pd

import furious_fastas as ff
from pep2prot.graph_ops import (get_minimal_protein_group_coverage_fast,
                                make_protein_peptide_graph)

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

peptides = list(set(diann_report_extract.peptide))
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
