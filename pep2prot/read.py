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


#TODO: modify furious_fastas to result in a quicker read in.
#TODO: use it aboves
def read_fastas(path, observed_prots):
    """Read in fastas files.

    Args:
        path (Path or str): Path to the fasta file used to sequence MS signals.
        observed_prots (set): Found proteins.
    Returns:
        dict: Mapping from proteins to their sequences.
    """
    fastas = UniprotFastas()
    fastas.read(path)
    fastas_df = pd.DataFrame.from_records(
        (' '.join(f.header.split(' ')[2:]),
         str(f),
         f.header.split(' ')[1]) 
        for f in fastas)
    fastas_df.columns = ['description', 'prot_seq', 'prot']
    fastas_df = fastas_df.set_index('prot')
    fastas_df['seq_len'] = fastas_df.prot_seq.map(len)
    assert all(fastas_df.groupby('prot').size() == 1), "The provided accesions are not unique and cannot form an index."
    assert all(r in fastas_df.index for r in observed_prots), "It seems that you are using a different set of fastas than the peptide annotation software before. Repent please."
    return fastas_df
