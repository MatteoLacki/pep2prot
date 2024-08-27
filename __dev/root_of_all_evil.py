%load_ext autoreload
%autoreload 2
from collections import Counter

import pandas as pd
from numba_progress import ProgressBar
from tqdm import tqdm

import ahocorasick
import furious_fastas as ff
from pep2prot.graph_ops import (get_adjacency_matrix,
                                get_minimal_protein_group_coverage)
from pep2prot.readers import read_DiaNN

fastas_path = "data/betterwheat/20211231_UniProtKB2021_04_4565TriticumAestivum_143648entries_172contaminants.fasta"

fastas = [(f.header.split(" ", 1)[0][1:], f.sequence) for f in ff.fastas(fastas_path)]
diann_report_extract = pd.read_parquet("data/diann_report_extract.parquet")

covering_protein_groups, adjacency_matrix = get_minimal_protein_group_coverage(
    peptides=list(set(diann_report_extract.peptide)),
    fastas=fastas,
    min_number_of_peptides=3,
)


protein_sequences = [sequence for header, sequence in fastas]
big_fasta = " ".join(protein_sequences)
peptides = list(set(diann_report_extract.peptide))


peps_found_mask = adjacency_matrix.any(axis=0)
pep_ids_found = set(map(int,peps_found_mask.nonzero()[0]))

automaton = ahocorasick.Automaton()
for pep_id,pep in enumerate(peptides):
    automaton.add_word(pep, pep_id)
automaton.make_automaton()
test2 = list(automaton.iter(big_fasta))

for end, pep_id in test2:
    pep = peptides[pep_id]
    start = end - len(pep) + 1
    assert big_fasta[start:end+1] == pep

pep_ids_found_cnts = Counter((pep_id for end, pep_id in test2))
assert set(pep_ids_found_cnts) == pep_ids_found, "Results of naive implementation and text match do not aggree."

# OK, but what we still need now is to construct a graph out of it all.

# so we need to index the fasta.
def is_nondecreasing(xx):
    xx = iter(xx)
    x_prev = next(xx)
    for x in xx:
        if x < x_prev:
            return False
        x_prev = x
    return True

is_nondecreasing(a for a,b in test2)
