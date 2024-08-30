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


def get_false_discovery_proportion_estimates(scores: npt.NDArray|pd.Series, decoys: npt.NDArray|pd.Series) -> npt.NDArray:
    """
    ASSUMPTION: HIGHER SCORE IS BETTER.
    """
    assert len(scores) == len(decoys), "Provide the same number of scores as their decoy labels."
    df = pd.DataFrame({"score": scores, "decoy": decoys}).reset_index()
    df = df.sort_values("score", ascending=False)

    df.score = df.score.astype(np.float32)
    df.decoy = df.decoy.astype(np.bool_)

    decoy_cumsum = df.decoy.cumsum()
    target_cumsum = (~df.decoy).cumsum()

    df["estim_FDP"] = np.minimum(1, (decoy_cumsum + 1)/ target_cumsum)
    df = df.sort_values("index")

    return df.estim_FDP.to_numpy()

psm_report["estim_FDP"] = get_false_discovery_proportion_estimates(psm_report.score, psm_report.decoy)

scores, decoys = [0,0,0],[1,1,1]
get_false_discovery_proportion_estimates([0,0,0],[0,1,0])
get_false_discovery_proportion_estimates([0,0,0],[1,1,0])
get_false_discovery_proportion_estimates([0,0,0],[0,0,0])




psm_report = psm_report.sort_values("score")

def get_q_values(fdps):
    min_fdp = fdps[0]
    for fdp in fdps:
        min_fdp = min(min_fdp, fdp)
        yield min_fdp

psm_report["qvalue"] = list(get_q_values(psm_report.FDP))
psm_report = psm_report.reset_index().sort_values("index").set_index("index")

psm_report.sort_values("qvalue")
