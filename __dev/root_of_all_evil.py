%load_ext autoreload
%autoreload 2
import pandas as pd
from numba_progress import ProgressBar

import ahocorasick
import furious_fastas as ff
from pep2prot.graph_ops import (get_adjacency_matrix,
                                get_minimal_protein_group_coverage)
from pep2prot.readers import read_DiaNN

fastas_path = "data/betterwheat/20211231_UniProtKB2021_04_4565TriticumAestivum_143648entries_172contaminants.fasta"

fastas = [(f.header.split(" ", 1)[0][1:], f.sequence) for f in ff.fastas(fastas_path)]
diann_report_extract = pd.read_parquet("data/diann_report_extract.parquet")

covering_protein_groups = get_minimal_protein_group_coverage(
    peptides=list(set(diann_report_extract.peptide)),
    fastas=fastas,
    min_number_of_peptides=3,
)

protein_sequences = [sequence for header, sequence in fastas]

automaton = ahocorasick.Automaton()
big_fasta = "_".join(protein_sequences)
