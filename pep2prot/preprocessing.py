import pandas as pd
import numpy as np


def simplify_mods(mods):
    if mods and not mods is np.nan:
        mods = mods.replace(', ',',')
        mods_d = {}
        for mod in mods.split(','):
            m, pos = mod.split(' ')
            pos = pos[0] + pos[2:-1]
            prev_pos = mods_d.get(m,'')
            mods_d[m] = prev_pos+"_"+pos if prev_pos else pos
        return " ".join(m+"_"+pos for m, pos in mods_d.items())
    else:
        return ''

simplify_mods = np.vectorize(simplify_mods)

def trivial_mods_simplification(mods_column):
    return mods_column.str.replace(' C\\(', '_C(').str.replace('\\,','')

def preprocess_isoquant_peptide_report(D, mods_simplifier=simplify_mods):
    I_cols = [c for c in D.columns if "intensity in" in c]
    D[I_cols] = D[I_cols].fillna(0)
    D.rename(columns={'pre_homology_entries':'prots_str'}, inplace=True)
    D.modifier = mods_simplifier(D.modifier)
    D['pep'] = np.where(D.modifier, D.sequence + " " + D.modifier, D.sequence)
    D['prots'] = D.prots_str.str.split(',').apply(frozenset)
    D_pep = D.groupby('pep')
    assert np.all(D_pep.prots.nunique() == 1), "Different proteins appear to explain the same peptides in different clusters. How come? Repent."    
    D.set_index('pep', inplace=True)
    return D, I_cols


def complex_cluster_buster(D, I_cols, unique_columns, max_rt_deviation=1):
    """Merge the same peptides that were in different clusters.

    Filter out peptides too far away in retention time from the top scoring cluster.

    Args:
        D (pd.DataFrame): Indexed by pep.
        I_cols (iterable): Names of columns containing peptide intensities.
        unique_columns (iterable): Names of columns with peptide-specific features.
    Returns:
        pd.DataFrame : filtered data uniquely indexed by peptides.
    """
    DD = D.copy()
    D_pep = DD.groupby('pep')
    pep_size = D_pep.size()
    D_mul_uni = DD.groupby(pd.Series(np.where(pep_size.values > 1, 'mul', 'uni'), pep_size.index))
    uni = D_mul_uni.get_group('uni')
    mul = D_mul_uni.get_group('mul').copy()
    mul.sort_index(inplace=True)
    mul_pep = mul.groupby(mul.index)
    mul_scores = mul.peptide_annotated_max_score - mul_pep.peptide_annotated_max_score.max()
    mul_rt2tophit = mul.signal_rt - mul.signal_rt[mul_scores == 0]
    mul = mul[np.abs(mul_rt2tophit) < max_rt_deviation]
    mul_pep = mul.groupby(mul.index)
    mul_intensities = mul_pep[I_cols].sum()
    mul_descriptors = mul_pep[unique_columns].head(1)
    mul_res = pd.concat([mul_intensities, mul_descriptors], axis=1, sort=True)    
    res = pd.concat([uni[I_cols + unique_columns], mul_res], sort=True)
    return res


def simple_cluster_buster(D, I_cols, unique_columns):
    """Merge the same peptides that were in different clusters.

    The aggregation does not take into account any filtering and is dead simple.

    Args:
        D (pd.DataFrame): indexed by pep.
        I_cols (iterable): Names of columns containing peptide intensities.
        unique_columns (iterable): Names of columns with peptide-specific features.
    Returns:
        pd.DataFrame : filtered data uniquely indexed by peptides.
    """
    D_pep = D.groupby('pep')
    aggregated_intensities = D_pep[I_cols].sum()
    no_change_here = D_pep[unique_columns].head(1)
    res = pd.concat([aggregated_intensities, no_change_here], axis=1, sort=True)
    return res

