import multiprocessing as mp

import furious_fastas as ff
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pep2prot.graph_ops import (
    check_protein_groups_are_OK,
    get_adjacency_matrix,
    get_protein_group_cover,
    get_protein_groups,
    sparsify_adjacency_martrix,
)
from pep2prot.readers import SimpleReplace

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
assert check_protein_groups_are_OK(
    protein_groups, adjacency_matrix
), "Some proteins in the same protein group do not share the same peptides (pattern of 0s and 1s differs)."

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
protein_group_cover = get_protein_group_cover(prot_pep_edges, cpu_cnt=cpu_cnt)

# comment: of course, merging peptides makes absolutely no sense in terms of a set cover. So we are right not to do so.
# what now?
# * Missing some form of reporting: for each peptide get a
# * More adapters for different peptide reports
# * Add the q-value calculations as preprocessing:
#   * we need ion level ones
# *
