import pandas as pd

from aa2atom import aa2atom, atom2mass
from aa2atom.aa2atom import UnknownAminoAcid


def get_stats(prots_min_I, prots_I, prots_max_I):
    """Get some stats on the deconvolution process.

    Args:
        prots_min_I (pd.DataFrame): The minimal assigned protein intensities.
        prots_I (pd.DataFrame): The deconvoluted protein intensities.
        prots_max_I (pd.DataFrame): The maximal assigned protein intensities.
    Returns:
        pd.DataFrame: a data frame summarizing the frequencies of particular conditions enumerated in the columns.
    """
    res = pd.concat([(prots_min_I < prots_I).sum()/len(prots_I),
                     (prots_I < prots_max_I).sum()/len(prots_I),
                     (prots_min_I == prots_I).sum()/len(prots_I),
                     (prots_I == prots_max_I).sum()/len(prots_I)],
                    axis=1)
    res.columns=['min < dec', 'dec < max', 'min = dec', 'dec = max']
    return res



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



def get_info_on_prots(prot_groups, fastas):
    """Get more information on the protein groups.

    Get the representative protein, its description, sequence, mass, etc.
    For a given group of proteins, the protein with the shortest sequence is chosen. In case of draw, the least massive protein is chosen. In case of a second draw (e.g. in case of isoforms), the protein with name earlier in the lexicographic order is chosen.

    Args:
        prot_groups (iterable): Consecutive diffrent protein groups.
        fastas (pd.Dataframe): Info on fastas.
    Returns:
        pd.DataFrame: mapping between protein group and its representative.
    """
    trivial_prot_reps = pd.DataFrame((rg,r) for rg in prot_groups
                                     if len(rg)==1 for r in rg)
    trivial_prot_reps.columns = ('prot', 'representative protein')
    trivial_prot_reps['protein_group'] = trivial_prot_reps['representative protein']
    trivial_prot_reps = trivial_prot_reps.join(fastas,
                                               on='representative protein')
    trivial_prot_reps['monoisotopic mass'] = trivial_prot_reps.sequence.map(aa2mass)

    intriguing_prot_reps = pd.DataFrame((rg,r) for rg in prot_groups
                                        if len(rg) > 1 for r in rg)
    intriguing_prot_reps.columns = ('prot', 'representative protein')
    intriguing_prot_reps = intriguing_prot_reps.join(fastas,
                                                     on='representative protein')
    intriguing_prot_reps['monoisotopic mass'] = intriguing_prot_reps.sequence.map(aa2mass)
    intriguing_prot_reps.sort_values(['prot','seq_len','monoisotopic mass','representative protein'],
                                     inplace=True)
    intriguing_prot_reps = intriguing_prot_reps.groupby('prot').head(1)
    intriguing_prot_reps['protein_group'] = intriguing_prot_reps.prot.map(lambda x: " ".join(x))
    res = pd.concat([trivial_prot_reps,
                     intriguing_prot_reps[trivial_prot_reps.columns]])
    res.set_index('prot', inplace=True)
    return res
