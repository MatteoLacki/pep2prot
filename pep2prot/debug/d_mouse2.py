%load_ext autoreload
%autoreload 2

from pathlib import Path
import pandas as pd

from pep2prot.analyse import isoquant_peptide_report
from pep2prot.read import read_isoquant_peptide_report
from pep2prot.preprocessing import complex_cluster_buster

test_data = Path("/Users/matteo/Projects/pep2prot/pep2prot/data/mouse")
pep_rep_path = test_data/''
fastas_path = test_data/'mouse.fasta'
