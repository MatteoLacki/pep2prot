import pandas as pd
import numpy as np


def preprocess_isoquant_peptide_report(D):
    D.modifier = D.modifier.fillna('')
    I_cols = [c for c in D.columns if "intensity in" in c]
    D[I_cols] = D[I_cols].fillna(0)
    D.rename(columns={'pre_homology_entries':'prots'}, inplace=True)
    D['pep'] = np.where(D['modifier'], D['sequence'] + "_" + D['modifier'], D['sequence'])
    D.prots = D.prots.str.split(',').apply(frozenset)
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
    D_pep = DD.groupby(DD.index)
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
    D_pep = D.groupby(D.index)
    aggregated_intensities = D_pep[I_cols].sum()
    no_change_here = D_pep[unique_columns].head(1)
    res = pd.concat([aggregated_intensities, no_change_here], axis=1, sort=True)
    return res

