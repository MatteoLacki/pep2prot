from pathlib import Path
import pandas as pd

from furious_fastas.fastas import UniprotFastas


def read_isoquant_peptide_report(path):
    """Read an IsoQuant peptide report.

    Args:
        path (Path or str): path to the file with the report.
    Return:
        pd.DataFrame: raw IsoQuant peptide report.
    """
    path = Path(path).expanduser()
    return pd.read_csv(path, encoding = "ISO-8859-1")


def read_n_check_fastas(path, observed_prots):
    """Read in fastas files.

    Args:
        path (Path or str): Path to the fasta file used to sequence MS signals.
        observed_prots (set): Found proteins.
    Returns:
        dict: Mapping from proteins to their sequences.
    """
    fastas = UniprotFastas()
    fastas.read(path)
    prot2seq = {f.header.split(' ')[1]: str(f) for f in fastas}
    assert all(r in prot2seq for r in observed_prots), "It seems that you are using a different set of fastas than the peptide annotation software before. Repent please."
    return {r:f for r,f in prot2seq.items() if r in observed_prots}