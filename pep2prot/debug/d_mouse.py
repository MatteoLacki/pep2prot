%load_ext autoreload
%autoreload 2

from pathlib import Path
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 30)#display whole column without truncation
from pandas import DataFrame as df
import numpy as np

from pep2prot.analyse import isoquant_peptide_report
from pep2prot.read import read_isoquant_peptide_report
from pep2prot.preprocessing import complex_cluster_buster

test_data = Path(r"~/Projects/pep2prot/tests/sabine4").expanduser()
pep_rep_path = test_data/'2019-023 Wolf MAM fraction_user designed 20190425-161116_peptide_quantification_report.csv'
fastas_path = test_data/'mouse.fasta'



prot_I_nice, all_prot_nice, G, H, lonely, unsupported, rejected, prot_min_I, prot_I, prot_max_I, D, W, pep_I = isoquant_peptide_report(pep_rep_path, fastas_path, full_outcome=True)


r = 'CAC1A_MOUSE'
R = next(R for R in H.prots() if r in R)
prot_I_nice.loc[r]
D[[r in R for R in D.prots]]
H[R]
D.loc['TALDIK']
wrong_I = ['intensity in 2019-023-13 Wolf S3 1', 'intensity in 2019-023-14 Wolf S4 1',
           'intensity in 2019-023-15 Wolf S7 1', 'intensity in 2019-023-16 Wolf S8 1']
TALDIK_I = D.loc['TALDIK', wrong_I]
prot_I.loc[[R], wrong_I].iloc[0]

H.component(R).draw(with_labels=True, font_size=6)
P = frozenset({'YWASLR', 'TMALYNPIPVR', 'TALDIK'})

prot_around_R = list({R2 for P in H[R] for R2 in H[P]})
prot_I.loc[prot_around_R, wrong_I]
W[wrong_I].loc[prot_around_R]

np.finfo(type(W.iloc[0,0])).eps# machine precision

## More info needed!
pep_I
H

# pandas gets confused, if the objects used for ids are frozensets, i.e. iterable objects.
prot_maxI.loc[[R], wrong_I] # the max intensities are correct!
prot_minI.loc[[R], wrong_I]

Pdeg1 = [P for P in H [R] if H.deg(P)==1][0]
D.loc[Pdeg1, wrong_I] 