%load_ext autoreload
%autoreload 2

import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 30)#display whole column without truncation
from pandas import DataFrame as df
from pathlib import Path
import platform

from pep2prot.analyse import isoquant_peptide_report
from pep2prot.preprocessing import complex_cluster_buster



test_data = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()
pep_rep_path = test_data/'hye_peprep.csv'
fastas_path = test_data/'HYE.fasta'
cluster_buster = complex_cluster_buster

prot_I_nice, all_prot_nice, G, H, lonely, unsupported, rejected, prot_min_I, prot_I, prot_max_I, D, W, pep_I = \
 	isoquant_peptide_report(pep_rep_path, fastas_path, verbose=True, full_outcome=True)

D.index.name = 'peptide sequence'
I_cols = [c for c in D.columns if "intensity in " in c]
E = D[['prots'] + I_cols]
E.columns = [c.replace('intensity in 2019-019-0','I').replace('4','').replace('5','').replace('prots','protein') for c in E.columns]

pd.set_option('display.expand_frame_repr', True)
pd.set_option('display.max_colwidth', 50)

F = prot_I_nice[I_cols]
F.columns = E.columns[1:]

F.loc['IF5A1_HUMAN'][0:3].mean()
F.loc['IF5A1_HUMAN'][3:6].mean()


# T = \frac{\bar X - \bar Y}{ S \sqrt{(\frac{1}{n_X} - \frac{1}{n_Y}) }}
# S = \sqrt{\frac{(n_X - 1)s_X^2 + (n_Y - 1)s_Y^2}{n_X + n_Y - 2}}
p = r"/Volumes/GREEN/2019-019_HYE/2019-019 HYE neu_user designed 20190621-083010_quantification_report_Jorg.xlsx"
J = pd.read_excel(p, skiprows=1)
J = J.set_index('entry')
J = J[[c.replace('intensity in ','') for c in I_cols]]
J.isna().sum(axis=0)

