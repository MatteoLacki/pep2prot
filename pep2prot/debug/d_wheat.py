%load_ext autoreload
%autoreload 2

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotnine import *
from platform import system
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 30)#display whole column without truncation
# pd.set_option('display.expand_frame_repr', True)

if system() == 'Linux':
    p = r"/home/matteo/Projects/pep2prot/data/20190719_wheat/"
else:
    p = r"/Volumes/GREEN/20190719_wheat-files-for-testing-protein-inference-by-ML/"

p += r"2019-015_triticum-DB_user designed 20190713-160429_peptide_quantification_report.csv"
M = pd.read_csv(p, encoding = "ISO-8859-1")
M = M[[c for c in M.columns if 'intensity in ' in c]]

y = M.isna().sum(axis=1)
x = M.median(axis=1)
plt.scatter(y, np.log(x))

D = pd.DataFrame({'NA_cnt': y, 'I':x})
D.NA_cnt = pd.Categorical(D.NA_cnt, ordered=True)
(ggplot(D) +
 geom_point(aes(y='NA_cnt', x='I')) + 
 scale_x_log10() +
 xlab('Intensity') + 
 ylab('Missing Observations') + 
 theme_minimal())


(ggplot(D) +
 geom_boxplot(aes(x='NA_cnt', y='I')) +
 scale_y_log10()+
 xlab('Missing Observations')+
 ylab('Intensity')+
 theme_minimal())
