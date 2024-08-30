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
from pep2prot.graph_ops2 import get_minimal_protein_group_coverage

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

%%time
(
    protein_groups,
    proteins,
    graph,
    pep_protgroup_graph,
    covers,
    peptide_view,
) = get_minimal_protein_group_coverage(
    peptides = list(set(diann_report_extract.peptide)),
    fastas=fastas,
    _debug=True
)

len(set(diann_report_extract.PG))

diann_report_extract["PGset"] = diann_report_extract.PG.map(lambda pg: tuple(sorted(set(pg.split(";")))))

protein_groups["PGset"] = protein_groups.accessions.map(lambda pg: tuple(sorted([acc.split("|",3)[1] for acc in pg])))


len(set(diann_report_extract.PGset) - set(protein_groups.PGset))
len(set(protein_groups.PGset) - set(diann_report_extract.PGset))
len(set(protein_groups.PGset) & set(diann_report_extract.PGset))
