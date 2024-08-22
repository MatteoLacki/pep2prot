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


peptides = list(set(ms2rescore_input["raw_sequence"]))
proteins = pd.DataFrame({"sequence": [f.sequence for f in fastas]})
proteins.index.name = "prot_id"
adjacency_matrix = get_adjacency_matrix(list(proteins.sequence), peptides)
proteins["peptide_cnt"] = adjacency_matrix.sum(axis=1)
proteins["group"] = get_protein_groups(np.packbits(adjacency_matrix, axis=-1))
assert check_protein_groups_are_OK(
    proteins.group.to_numpy(),
    adjacency_matrix,
), "Some proteins in the same protein group do not share the same peptides (pattern of 0s and 1s differs)."

proteins_with_enough_peptides = proteins[proteins.peptide_cnt >= min_number_of_peptides]

# this is wrong! it should be that the groups should dictate what goes in.
# proteins.group[proteins.peptide_cnt >= min_number_of_peptides]

protein_groups_representatives = (
    proteins_with_enough_peptides.reset_index().groupby("group").first()
)
representitives_adjacency_matrix = adjacency_matrix[
    protein_groups_representatives.prot_id
]
assert np.all(
    protein_groups_representatives.peptide_cnt
    == representitives_adjacency_matrix.sum(axis=1)
), "Quality check failed."

prot_pep_edges = sparsify_adjacency_martrix(
    representitives_adjacency_matrix, protein_groups_representatives.prot_id.to_numpy()
)
protein_group_cover = get_protein_group_cover(prot_pep_edges, cpu_cnt=cpu_cnt)


# comment: of course, merging peptides makes absolutely no sense in terms of a set cover. So we are right not to do so.
# what now?
# * Missing some form of reporting: for each peptide get a
# * More adapters for different peptide reports
# * Add the q-value calculations as preprocessing:
#   * we need ion level ones
# *
