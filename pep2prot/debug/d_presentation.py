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
from plotnine import *

from pep2prot.read import read_isoquant_peptide_report, read_fastas
from pep2prot.preprocessing import preprocess_isoquant_peptide_report, get_protein_coverages, complex_cluster_buster, simple_cluster_buster 
from pep2prot.graphs.pep_prot_graph import ProtPepGraph
from pep2prot.intensities import get_prot_intensities
from pep2prot.df_ops import sum_top
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
assert DD.index.is_unique, "Some peptides are still across different clusters."

DD.index.name = 'peptide sequence'
EE = DD[['prots'] + I_cols]
EE.columns = [c.replace('intensity in 2019-019-0','I').replace('4','').replace('5','').replace('prots','protein') for c in EE.columns]

pd.set_option('display.expand_frame_repr', True)
pd.set_option('display.max_colwidth', 50)
EE