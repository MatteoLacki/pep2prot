%load_ext autoreload
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
observed_prots = {r for rg in D.prots for r in rg}
assert all(r in fastas.index for r in observed_prots), "It seems that you are using a different set of fastas than the peptide annotation software before. Repent please."
prots = get_protein_coverages(D, fastas)
uni_cols = ['peptide_overall_max_score','peptide_fdr_level',
        'peptide_overall_replication_rate','prots',
            'pre_homology_accessions','pi','mw']

DD = cluster_buster(D, I_cols, uni_cols) # agg same peptides in various clusters
assert np.all(DD.groupby(['pep','prots']).size() == 1), "Some peptides are still across different clusters."

G = ProtPepGraph((r,p) for p, rg in zip(DD.index, DD.prots) for r in rg)
lonely, unsupported = G.remove_lonely_and_unsupported(2)
H, rejected = G.get_minimal_graph()

pep2pepgr = {p:pg for pg in H.peps() for p in pg}
DDinH = DD.loc[pep2pepgr] # peps in H: no simple prot-pep pairs, no unnecessary prots?
DDinH['pepgr'] = DDinH.index.map(pep2pepgr)
pep_I = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities

prot_minI, prot_I, prot_maxI = get_prot_intensities(H, pep_I)
prot_info = summarize_prots(H, fastas, prots.pep_coverage)

plt.hist(prot_info.pep_coverage, bins=100)
plt.show()

prot_Inice = prettify_protein_informations(prot_I, prot_info)
prot_stats = get_stats(prot_minI, prot_I, prot_maxI)

prot_Inice.columns
any(';' in d for d in prot_Inice.description)

all_prots = get_full_report(prot_minI, prot_I, prot_maxI)
all_prots_nice = prettify_protein_informations(all_prots, prot_info)



