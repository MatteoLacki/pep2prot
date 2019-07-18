%load_ext autoreload
%autoreload 2

from pathlib import Path
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import networkx as nx
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 60)#display whole column without truncation

from aa2atom import aa2atom, atom2mass
from aa2atom.aa2atom import UnknownAminoAcid
from pep2prot.read import read_isoquant_peptide_report, read_n_check_fastas
from pep2prot.graphs import ProtPepGraph, BiGraph, get_peptide_protein_graph
from pep2prot.preprocessing import preprocess_isoquant_peptide_report, complex_cluster_buster, simple_cluster_buster

eps = np.finfo(type(weights.iloc[0,0])).eps# machine precision

min_pepNo_per_prot = 2
max_rt_deviation = 1
path = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()

D = read_isoquant_peptide_report(path/'peptide_report.csv')
D, I_cols = preprocess_isoquant_peptide_report(D)
# X = D.groupby(D.index).nunique().nunique() # the double unique beast!
unique_columns = ['peptide_overall_max_score','peptide_fdr_level','peptide_overall_replication_rate','prots','pre_homology_accessions','pi','mw']
D2 = complex_cluster_buster(D, I_cols, unique_columns, max_rt_deviation)
D3 = simple_cluster_buster(D, I_cols, unique_columns)
DD = D2

prot2seq = {r for rg in D2.prots for r in rg}
prot2seq = read_n_check_fastas(path/'mouse.fasta', prot2seq)

# aggregate the same peptides in different clusters
H2, RWEP2, BRR2 = get_peptide_protein_graph(D2)
H3, RWEP3, BRR3 = get_peptide_protein_graph(D3)
H = H2

# NOW: FINALLY THE BLOODY REPORT
pep2pepgr = {p:pg for pg in H.peps() for p in pg}
D2 = DD.loc[pep2pepgr] # peps in H (no simple prot-pep pairs)
D2.drop(('prots'), axis=1, inplace=True)
peps_I = D2[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities

prot_pep = ('prot','pep')
Hdf = pd.DataFrame.from_records(H.prot_pep_pairs(), columns=prot_pep, index=prot_pep)
peps2prots_max_I   = Hdf.join(peps_I, on='pep')
peps2prots         = peps2prots_max_I.index
##############################################################################
prots_max_I        = peps2prots_max_I.groupby(level='prot').sum()## RESULT ###
##############################################################################
unipeps            = pd.Index((pg for pg in H.peps() if H.degree(pg) == 1),
                              name='pep') # peps unique for some proteins
unipeps2prots_I    = peps2prots_max_I[peps2prots.get_level_values('pep').isin(unipeps)]
# there is one unique peptide per protein, so _I == _max_I == _min_I
uniprots_min_I     = unipeps2prots_I.reset_index('pep', drop=True)
uniprots           = uniprots_min_I.index # prots with some uniquely identified peptides
# blade = razor??? probably not, so stick with blade.
bladeprots         = pd.Index((r for r in H.prots() if not r in uniprots), name='prot')
bladeprots_zero_I  = pd.DataFrame(np.zeros(shape=(len(bladeprots),
                                                  len(unipeps2prots_I.columns))),
                                  index=bladeprots, columns=unipeps2prots_I.columns)
# prots = uniprots ⊔ bladeprots
#########################################################################
prots_min_I = pd.concat([uniprots_min_I, bladeprots_zero_I])## RESULT ###
#########################################################################
sorted_prots = np.sort(prots_min_I.index)
prots_min_I = prots_min_I.loc[sorted_prots]
prots_max_I = prots_max_I.loc[sorted_prots]
assert np.all(prots_min_I <= prots_max_I), "Some minimal intensities are not smaller then the maximal ones. Report to Matteo."

uniprots_curr_I = uniprots_min_I# reusable: intensity from unique peps pushed to uniprots
bladepeps = pd.Index({p for r in bladeprots 
                         for p in H[r] if all(rr in bladeprots for rr in H[p])},
                     name='pep')# peps with prots without unique peps
bladepeps_I = peps_I.loc[bladepeps]
bladepeps_protcnts = pd.Series((H.degree(p) for p in bladepeps), index=bladepeps)
# distribute blade peps intensity uniformly among neighbor prots: hence the division below
bladepeps2bladeprots_I = pd.DataFrame.from_records(
    ((r,p) for p in bladepeps for r in H[p]), 
    columns=prot_pep,
    index=prot_pep).join( bladepeps_I.div(bladepeps_protcnts, axis='index'), on='pep')
bladeprots_curr_I = bladepeps2bladeprots_I.groupby(level='prot').sum()
# bladeprots might receive some more intensity from other peptides, hence _curr_I.
# peps = unipeps ⊔ bladepeps ⊔ otherpeps
# unipeps neighbor uniprots, bladepeps neighbor bladeprots, otherpeps neighbor some prots from both sets
otherpeps    = pd.Index({p for p in H.peps() if p not in unipeps and p not in bladepeps}, name='pep')
otherpeps_I  = peps_I.loc[otherpeps]
# above intensities will be distributed proportionally to intensities prots received from bladepeps and unipeps.
# call these the current intensities, curr_I.
# bladeprots and uniprots are disjoint sets, so we concat them
prots_curr_I = pd.concat([uniprots_curr_I, bladeprots_curr_I])

# populating edges with other intensities from otherpeps. Let's call their prots mixprots
otherpeps2mixprots_I = pd.DataFrame(index=peps2prots).join(otherpeps_I, on='pep', how='right')
otherpeps2mixprots = otherpeps2mixprots_I.index
mixprots = otherpeps2mixprots_I.index.get_level_values('prot').unique()
# need weights for otherpeps2mixprots_I
weights = pd.DataFrame(index=otherpeps2mixprots).join(prots_curr_I, on='prot')
weights += eps
weights_mixprot_I = weights.groupby(level='pep').sum()
weights = weights.div(weights_mixprot_I, axis='index')

assert not np.any(weights.isnull()), "Weight cannot result in any NaNs."
otherpeps2mixprots_I = otherpeps2mixprots_I * weights
#update only mixprots
prots_curr_I.loc[mixprots] += otherpeps2mixprots_I.groupby(level='prot').sum()
########################################################################
prots_I = prots_curr_I.loc[sorted_prots] ## RESULT #####################
########################################################################
assert np.all(prots_min_I <= prots_I), "Some deconvoluted intensities are smaller then minimal intensities."
assert np.all(prots_I <= prots_max_I), "Some deconvoluted intensities are larger then maximal intensities."

# some stats needed on equalities
def get_stats(prots_min_I, prots_I, prots_max_I):
    res = pd.concat([(prots_min_I < prots_I).sum()/len(prots_I),
                      (prots_I < prots_max_I).sum()/len(prots_I),
                      (prots_min_I == prots_I).sum()/len(prots_I),
                      (prots_I == prots_max_I).sum()/len(prots_I)],
                    axis=1)
    res.columns=['min < dec', 'dec < max', 'min = dec', 'dec = max']
    return res





# reporting finally
for rgr in R.prots():
    if len(rgr) > 2:
        break

def aa2mass(aa, which_mass='monoisotopic', _big_error_mass=1e12):
    try:
        f = aa2atom(aa)
        return atom2mass(f, which_mass)
    except UnknownAminoAcid as uaa:
        return _big_error_mass

# maybe these operations could be performed beforand on all things?
# select all the proteins in R and build up one dataframe.

# not all needed: only those, that ain't unique: and maybe not those, that are within

r_df = pd.DataFrame.from_records(((r,i,c) for i,c in enumerate(R.components())
                                           for r_gr in c.prots() if len(r_gr)>1
                                            for r in r_gr),
                                 columns=('prot','cc_index','cc'))
# r_df = pd.DataFrame({'prot':list(r for r_gr in R.prots() if len(r_gr) > 1 for r in r_gr)})
r_df['seq'] = r_df.prot.map(accession2seq)
r_df['len'] = r_df.seq.map(len)
r_df['mass'] = r_df.seq.map(aa2mass) # this take the most runtime
r_df.sort_values(['cc_index','len', 'mass', 'prot'], inplace=True)



#TODO: rewrite without pandas and faster.
def trivial_representative_prot(r_gr):
    """Choose a representative of a group of proteins."""
    r_df = pd.DataFrame({'prot':list(r_gr)})
    r_df['seq'] = r_df.prot.map(accession2seq)
    r_df['seq_len'] = r_df.seq.map(len)
    r_df['mass'] = r_df.seq.map(aa2mass)
    r_df.sort_values(['seq_len', 'mass', 'prot'], ascending=[True, True, True], inplace=True)
    return dict(r_df.iloc[0])


X = [trivial_representative_prot(r_gr)['prot'] if len(r_gr) > 2 else next(r for r in r_gr)
     for r_gr in R.prots()]




it = R.components()
cc = next(it)
r = next(cc.B())
D = D.set_index('pep')

# fastes way to divide the data
pep2cc = {p: i for i, cc in enumerate(R.components()) for pg in cc.peps() for p in pg}
x = list(D.groupby(pep2cc))
x[0]

# add in the fasta file


def report(R):
    for cc in R.components():
        p_cnt, r_cnt = cc.nodes_cnt()
        if r_cnt == 0: # a single peptide
            yield 'gowno' 
        if r_cnt == 1:
            yield simple_report(cc)

