import pandas as pd
import numpy as np

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
                     (prots_I == prots_max_I).sum()/len(prots_I),
                     (prots_min_I == prots_max_I).sum()/len(prots_I)],
                    axis=1)
    res.columns=['min < dec', 'dec < max', 'min = dec', 'dec = max', 'min = max']
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


def summarize_prots(H, fastas, prot_coverages):
    """Get more information on the protein groups.

    Get the representative protein, its description, sequence, mass, etc.
    For a given group of proteins, the protein with the shortest sequence is chosen. In case of draw, the least massive protein is chosen. In case of a second draw (e.g. in case of isoforms), the protein with name earlier in the lexicographic order is chosen.

    Args:
        H (PepProtGraph): The estimated graph.
        fastas (pd.Dataframe): Info on fastas.
        prot_coverages (pd.Series or pd.DataFrame): Coverages for a superset of found proteins.
    Returns:
        pd.DataFrame: mapping between protein group and its representative.
    """
    found_prot_groups = set(H.prots())
    trivial_prot_reps = pd.DataFrame((rg,r) for rg in found_prot_groups
                                     if len(rg)==1 for r in rg)
    trivial_prot_reps.columns = ('prot', 'representative protein')
    trivial_prot_reps['protein_group'] = ''
    trivial_prot_reps = trivial_prot_reps.join(fastas,on='representative protein')
    trivial_prot_reps['monoisotopic mass'] = trivial_prot_reps.prot_seq.map(aa2mass)

    intriguing_prot_reps = pd.DataFrame((rg,r) for rg in found_prot_groups
                                        if len(rg) > 1 for r in rg)
    intriguing_prot_reps.columns = ('prot', 'representative protein')
    intriguing_prot_reps = intriguing_prot_reps.join(fastas,
                                                     on='representative protein')
    intriguing_prot_reps['monoisotopic mass'] = intriguing_prot_reps.prot_seq.map(aa2mass)
    intriguing_prot_reps.sort_values(['prot','seq_len','monoisotopic mass','representative protein'],
                                     inplace=True)
    intriguing_prot_reps = intriguing_prot_reps.groupby('prot').head(1)
    prot_groups = []
    for r,rg in zip(intriguing_prot_reps['representative protein'],
                    intriguing_prot_reps.prot):    
        set(rg).remove(r)
        prot_groups.append(" ".join(rg))
    intriguing_prot_reps['protein_group'] = prot_groups
    res = pd.concat([trivial_prot_reps,
                     intriguing_prot_reps[trivial_prot_reps.columns]])
    res = res.set_index('prot')
    res = res.join(prot_coverages, on='representative protein')

    prot2pep = pd.DataFrame.from_records(H.prot_pep_pairs(), columns=('prot','pep'))
    prot2pep['pep_cnt'] = prot2pep.pep.map(len)
    res['peptides count'] = prot2pep.groupby('prot').pep_cnt.sum()
    return res


def prettify_protein_informations(prot_intensities, prot_info):
    """Take intensities and make them more readable.

    Args:
        prot_intensities (pd.DataFrame): Indexed uniquely by proteins, columns contain intensities for each run.
        prot_info (pd.DataFrame): Indexed by precisely the same proteins as 'prot_intensities', columns contain information about the proteins.
    """
    X = prot_info.join(prot_intensities.applymap(int))
    X = X.set_index('representative protein')
    X.pep_coverage = ["{:.3%}".format(pc) for pc in X.pep_coverage]
    X = X.rename(columns={'prot_seq':'sequence', 'seq_len':'sequence length', 'pep_cover': 'peptide coverage', 'protein_group': 'other proteins in group'})
    X = X.drop(columns='sequence')
    X['monoisotopic mass'] = [round(m,3) for m in X['monoisotopic mass']]
    return X


def get_full_report(prots_min_I, prots_I, prots_max_I):
    """Prepare a full report out of partial reports.

    Args:
        prots_min_I (pd.DataFrame): The minimal estimated intensities.
        prots_I (pd.DataFrame): The deconvolved intensities.
        prots_max_I (pd.DataFrame): The maximal estimated intensities.

    Returns:
        pd.DataFrame: A DAIMOVIK of three different reports.
    """
    all_prots = np.zeros(dtype=prots_I.values.dtype, shape=(prots_I.shape[0],prots_I.shape[1]*3))
    all_prots[:,0::3] = prots_min_I.values
    all_prots[:,1::3] = prots_I.values
    all_prots[:,2::3] = prots_max_I.values
    columns = []
    for I in prots_min_I.columns:
        columns.append('minimal '+I)
        columns.append('deconvoled '+I)
        columns.append('maximal '+I)    
    return pd.DataFrame(all_prots, columns=columns, index=prots_I.index)