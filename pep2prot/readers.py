import re

import pandas as pd

import duckdb


class SimpleReplace:
    def __init__(self, pattern: str = r"\[.*?\]"):
        self._pattern = re.compile(pattern)

    def apply(self, string: str, withwhat: str = ""):
        return re.sub(self._pattern, withwhat, string)


def read_ms2rescore_peptide_report(path: str) -> pd.DataFrame:
    """
    Extract a peptide report from ms2rescore mokapot.peptides.txt input.
    """
    ms2rescore_input = pd.read_csv(path, sep="\t")
    ms2rescore_input = ms2rescore_input[
        [
            "peptide",
            "expmass",
            "retention_time",
            "charge",
            "mokapot q-value",
            "protein_list",
        ]
    ]
    mods_bracket_anihilator = SimpleReplace()
    ms2rescore_input["raw_sequence"] = ms2rescore_input.peptide.map(
        mods_bracket_anihilator.apply
    )
    return ms2rescore_input


def read_DiaNN(path: str) -> pd.DataFrame:
    """Read columns from the DiaNN 'report.tsv'."""
    return (
        duckdb.connect()
        .query(
            """
    SELECT DISTINCT 
    "Stripped.Sequence" AS peptide,
    "Protein.Group" AS PG,
    FROM '{diann_report}'
    """.format(
                diann_report=path
            )
        )
        .df()
    )
