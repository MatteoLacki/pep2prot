import pandas as pd

from aa2atom import aa2atom, atom2mass
from aa2atom.aa2atom import UnknownAminoAcid


#TODO: add this to the furious fastas module.
#TODO: this might not be the best idea to hide an exception as a very big mass.
def aa2mass(aa, which_mass='monoisotopic', _big_error_mass=1e12):
    """Calculate the mass of an aminoacidic sequence.

    Args:
        which_mass ('monoisotopic' or 'average'): Monoisotopic or average mass.
        _big_error_mass (float): Mass reported in case of an error.
    """
    try:
        f = aa2atom(aa)
        return atom2mass(f, which_mass)
    except UnknownAminoAcid as uaa:
        return _big_error_mass


def choose_reps(prot_groups, prot2seq):
    """Choose representatives of protein groups.

    For a given group of proteins, the protein with the shortest sequence is chosen. In case of draw, the least massive protein is chosen. In case of a second draw (e.g. in case of isoforms), the protein with name earlier in the lexicographic order is chosen.

    Args:
        prot_groups (iterable): Consecutive diffrent protein groups.
        prot2seq (dict): Maps proteins to their sequences.
    Returns:
        pd.DataFrame: mapping between protein group and its representative.
    """
    trivial_prot_reps = pd.DataFrame((rg,r) for rg in prot_groups if len(rg)==1 for r in rg)
    trivial_prot_reps.columns = ('protgroup', 'prot')
    trivial_prot_reps['protein_group'] = trivial_prot_reps.prot
    intriguing_prot_reps = pd.DataFrame((rg,r) for rg in prot_groups if len(rg) > 1 for r in rg)
    intriguing_prot_reps.columns = ('protgroup', 'prot')
    intriguing_prot_reps['seq'] = intriguing_prot_reps.prot.map(prot2seq)
    intriguing_prot_reps['seq_len'] = intriguing_prot_reps.seq.map(len)
    intriguing_prot_reps['mass'] = intriguing_prot_reps.seq.map(aa2mass)
    intriguing_prot_reps.sort_values(['protgroup','seq_len','mass','prot'], inplace=True)
    intriguing_prot_reps = intriguing_prot_reps.groupby('protgroup').head(1)
    intriguing_prot_reps['protein_group'] = intriguing_prot_reps.protgroup.map(lambda x: " ".join(x))
    res = pd.concat([intriguing_prot_reps[['protgroup', 'prot', 'protein_group']], trivial_prot_reps])
    res.columns = ['prot','repr','protein_group']
    res.set_index('prot', inplace=True)
    return res
