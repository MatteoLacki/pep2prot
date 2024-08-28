%load_ext autoreload
%autoreload 2
import multiprocessing as mp
import typing
from collections import Counter, defaultdict

import networkx as nx
import numba
import numpy as np
import numpy.typing as npt
import pandas as pd
from numba import types
from numba.typed import Dict
from numba_progress import ProgressBar
from tqdm import tqdm

psm_report = pd.read_parquet("data/SAGE_PSM_REPORT.parquet")

psm_report = psm_report.sort_values("score", ascending=False)


decoy_cumsum = np.cumsum(psm_report.decoy)
target_cumsum = np.cumsum(~psm_report.decoy)

psm_report["FDP"] = np.minimum(1, (decoy_cumsum + 1)/ target_cumsum)
psm_report = psm_report.sort_values("score")

def get_q_values(fdps):
    min_fdp = fdps[0]
    for fdp in fdps:
        min_fdp = min(min_fdp, fdp)
        yield min_fdp

psm_report["qvalue"] = list(get_q_values(psm_report.FDP))
psm_report = psm_report.reset_index().sort_values("index").set_index("index")