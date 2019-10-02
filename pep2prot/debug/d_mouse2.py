%load_ext autoreload
%autoreload 2

from pathlib import Path
import pandas as pd
from platform import system
import matplotlib.pyplot as plt
from math import log
import numpy as np

from pep2prot.analyse import isoquant_peptide_report
from pep2prot.read import read_isoquant_peptide_report
from pep2prot.preprocessing import complex_cluster_buster

if system()=='Linux':
    test_data = Path("/home/matteo/Projects/pep2prot/pep2prot/data/")
else:
    test_data = Path("/Users/matteo/Projects/pep2prot/pep2prot/data/")

pep_rep_path = test_data/'peptide_report.csv'
fastas_path = test_data/'mouse.fasta'
M = pd.read_csv(pep_rep_path, encoding = "ISO-8859-1")
M = M[[c for c in M.columns if 'intensity in ' in c]]

D = pd.DataFrame({'NA_cnt': M.isna().sum(axis=1), 'I':M.median(axis=1)})
D.NA_cnt = pd.Categorical(D.NA_cnt, ordered=True)

# p = (ggplot(D) +
#  geom_boxplot(aes(x='NA_cnt', y='I')) +
#  scale_y_log10()+
#  xlab('Missing Observations')+
#  ylab('Intensity')+
#  theme_minimal())

dM = M.sub(M.median(axis=1), axis='index')
y = np.array([i for c in dM for i in dM[c] if i and not pd.isna(i)])
plt.hist(y, bins=10000)

logX = np.log(M)
dlogX = logX.sub(logX.median(axis=1), axis='index')
x = np.array([i for c in dlogX for i in dlogX[c] if i and not pd.isna(i)])
plt.hist(x, bins=1000)


logX = np.log(M)
dlogX_mean = logX.sub(logX.mean(axis=1), axis='index')
z = np.array([i for c in dlogX_mean for i in dlogX_mean[c] if i and not pd.isna(i)])
plt.hist(z, bins=1000)


