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
from pep2prot.read import read_fastas

from furious_fastas.fastas import UniprotFastas


# the furious fastas have to take into account
path2 = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()
fasta_path2 = path2/'mouse.fasta'
fastas2 = UniprotFastas()
fastas2.read(fasta_path2)
fastas2[0]

path = Path('/home/matteo/Projects/pep2prot/tests/sabine2/Mouse_and_contaminats.fasta')
fastas = UniprotFastas()
fastas.read(path)
f = fastas[0]
f

f.header.split(' ')[1]

fastas_df = pd.DataFrame.from_records(
    (' '.join(f.header.split(' ')[2:]),
     str(f),
     f.header.split(' ')[1]) 
    for f in fastas)
fastas_df.columns = ['description', 'prot_seq', 'prot']
fastas_df = fastas_df.set_index('prot')
fastas_df['seq_len'] = fastas_df.prot_seq.map(len)

fastas_df.index.nunique()

assert all(fastas_df.groupby('prot').size() == 1), "The provided accesions are not unique and cannot form an index."
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

# from pep2prot.analyse import isoquant_peptide_report

# path         = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()
# pep_rep_path = path/'peptide_report.csv'
# fasta_path   = path/'mouse.fasta'
# prots_I_nice, all_prots_nice = isoquant_peptide_report(pep_rep_path, fasta_path)
# prots_I_nice.to_csv(Path(r"~/Projects/pep2prot/protein_report.csv").expanduser())
# all_prots_nice.to_csv(Path(r"~/Projects/pep2prot/full_protein_report.csv").expanduser())

