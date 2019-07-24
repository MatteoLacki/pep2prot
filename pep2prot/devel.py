%load_ext autoreload
%autoreload 2

from collections import Counter
import matplotlib.pyplot as plt

import networkx as nx
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 30)#display whole column without truncation
from pathlib import Path

from pep2prot.analyse import isoquant_peptide_report

path         = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()
pep_rep_path = path/'peptide_report.csv'
fasta_path   = path/'mouse.fasta'

prots_I_nice, all_prots_nice = isoquant_peptide_report(pep_rep_path, fasta_path)

prots_I_nice.to_csv(Path(r"~/Projects/pep2prot/protein_report.csv").expanduser())
all_prots_nice.to_csv(Path(r"~/Projects/pep2prot/full_protein_report.csv").expanduser())

# DEBUG
# prots_I
# error_maker = 'MPRIP_MOUSE'
# sample = "intensity in 2019-023-22 Wolf 1"
# ergr = [rg for rg in prots_I.index for r in rg if error_maker in r][0]
# er_peps = [p for pg in H[ergr] for p in pg]

# 35534/2
# H[frozenset(['IEDLQR'])]
# H.subgraph(nx.node_connected_component(H, ergr))

# def peep():
#     for pg in H[ergr]:
#         yield pg
#         for rg in H[pg]:
#             yield rg
# MARKOV_BLANKET = H.subgraph(set(peep()))
# MARKOV_BLANKET.draw(with_labels=True)

# D.loc[er_peps,sample]
# the_other_protein = list(MARKOV_BLANKET.prots())[1]
# prots_I_nice.loc[error_maker,sample]
# prots_I_nice.loc[the_other_protein,sample]