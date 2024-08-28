%load_ext autoreload
%autoreload 2
import pandas as pd

import furious_fastas as ff
from pep2prot.graph_ops import get_minimal_protein_group_coverage_fast

pd.set_option('display.max_columns', None)
fastas_path = "data/betterwheat/20211231_UniProtKB2021_04_4565TriticumAestivum_143648entries_172contaminants.fasta"

fastas = [(f.header.split(" ", 1)[0][1:], f.sequence) for f in ff.fastas(fastas_path)]
diann_report_extract = pd.read_parquet("data/diann_report_extract.parquet")

%%time
protein_groups, proteins, graph, protein_representatives_covers = get_minimal_protein_group_coverage_fast(
    peptides = list(set(diann_report_extract.peptide)),
    fastas = fastas,
    _debug = True,
)
