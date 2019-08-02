load_ext autoreload
%autoreload 2

import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 30)#display whole column without truncation
from pandas import DataFrame as df
from pathlib import Path
from collections import Counter
import networkx as nx
import matplotlib.pyplot as plt

from pep2prot.read import read_isoquant_peptide_report, read_fastas
from pep2prot.preprocessing import preprocess_isoquant_peptide_report, get_protein_coverages, complex_cluster_buster, simple_cluster_buster 
from pep2prot.graphs.pep_prot_graph import ProtPepGraph
from pep2prot.intensities import get_prot_intensities
from pep2prot.postprocessing import summarize_prots, get_stats, prettify_protein_informations, get_full_report

test_data = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()
pep_rep_path = test_data/'hye_peprep.csv'
fastas_path = test_data/'HYE.fasta'
cluster_buster = complex_cluster_buster

D = read_isoquant_peptide_report(pep_rep_path)
D, I_cols = preprocess_isoquant_peptide_report(D)
fastas = read_fastas(fastas_path)

pep = D.index
prot = D.prots
pepseq = D.sequence


peptides, proteins, peptide_sequences = D.index,D.sequence,D.prots

from pep2prot.string_ops import find_indices3
from pep2prot.range_ops import covered_area
from functools import partial

D[['sequence', 'prots']].reindex()


X = pd.DataFrame(((p,ps,r) for p,ps,rg in zip(D.index,D.sequence,D.prots) for r in rg),
                 columns=('pep','pepseq','prot'))
X['protseq']= X.prot.map(fastas.prot_seq)
X['sub_idx']= [find_indices3(p,psp) for p,psp in zip(X.pepseq, X.protseq)]
X.sub_idx   = X.sub_idx.map(frozenset)

sub_indices = X[['prot','sub_idx']].set_index('prot')
sub_indices = sub_indices.groupby('prot').sub_idx.apply(lambda x: sorted(set([]).union(*x)))
sub_indices.map(covered_area)

sorted_intervals = sub_indices[0]
covered_area(sorted_intervals)




prots = pd.DataFrame({'prot':X.prot.unique()}).set_index('prot').join(fastas)
prots['pep_cover_len'] = sub_indices.map(reverse_sorted_range_list_len)

sub_indices.iloc[0]
prots['pep_cover_len'][0]




assert np.all(prots.seq_len >= prots.pep_cover_len), "Peptide coverage exceeded protein length. If you ask me, if that is not an error, then what is?"

prots['pep_coverage'] = prots.pep_cover_len/prots.seq_len


