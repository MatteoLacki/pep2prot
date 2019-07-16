%load_ext autoreload
%autoreload 27

import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import networkx as nx
import numpy as np
from pathlib import Path
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 40)#display whole column without truncation

from furious_fastas.fastas import UniprotFastas
from furious_fastas.contaminants import uniprot_contaminants
from aa2atom import aa2atom, atom2mass
from aa2atom.aa2atom import UnknownAminoAcid
from pep2prot.graphs import ProtPepGraph, BiGraph

min_pepNo_per_prot = 2
max_rt_deviation = 1
path = r"~/Projects/pep2prot/pep2prot/data"

path = Path(path).expanduser()
D = pd.read_csv(path/'peptide_report.csv', encoding = "ISO-8859-1")
fastas = UniprotFastas()
fastas.read(path/'mouse.fasta')

# preprocessing
D.modifier = D.modifier.fillna('')
I_cols = [c for c in D.columns if "intensity in" in c]
D[I_cols] = D[I_cols].fillna(0)
D.rename(columns={'pre_homology_entries':'prots'}, inplace=True)
D['pep'] = np.where(D['modifier'], D['sequence'] + "_" + D['modifier'], D['sequence'])
D.prots = D.prots.str.split(',').apply(frozenset)
accession2seq = {f.header.split(' ')[1]: str(f) for f in fastas}
obs_accessions = {r for rg in D.prots for r in rg}
assert obs_accessions <= set(accession2seq), "It seems that you are using a different set of fastas than the peptide annotation software before. Repent please."
D_pep = D.groupby('pep')
assert np.all(D_pep.prots.nunique() == 1), "Different proteins appear to explain the same peptides in different clusters. How come? Repent."
# aggregate the same peptides in different clusters

# X = D_pep.nunique().nunique() # the double unique beast!

unique_columns = ['peptide_overall_max_score','peptide_fdr_level','peptide_overall_replication_rate','prots','pre_homology_accessions', 'pi', 'mw']


def complex_cluster_buster(D, I_cols, unique_columns):
    """Merge the same peptides that were in different clusters.

    Filter out peptides too far away in retention time from the top scoring cluster.

    Args:
        D (pd.DataFrame): indexed by pep.
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
    res.reset_index(inplace=True)
    return res
    

def simple_cluster_buster(D, I_cols, unique_columns):
    """Merge the same peptides that were in different clusters.

    The aggregation does not take into account any filtering and is dead simple.
    """
    D_pep = D.groupby(D.index)
    aggregated_intensities = D_pep[I_cols].sum()
    no_change_here = D_pep[unique_columns].head(1)
    res = pd.concat([aggregated_intensities, no_change_here], axis=1, sort=True)
    res.reset_index(inplace=True)
    res.rename(columns={'index':'pep'}, inplace=True)
    return res

E = complex_cluster_buster(D, I_cols, unique_columns)
F = simple_cluster_buster(D, I_cols, unique_columns)


# Fun begins here:   peptide BLUE   protein: RED
def get_graph(data, min_pepNo_per_prot=2):
    """Get the petide-protein-groups graph."""
    G = ProtPepGraph((r,p) for rs, p in zip(data.prots, data.pep) for r in rs)
    prots_without_enough_peps = [r for r in G.prots() if G.degree(r) < min_pepNo_per_prot]
    G.remove_nodes_from(prots_without_enough_peps)
    H = G.form_groups()
    HMC = H.greedy_minimal_cover() # Her Majesty's Minimal Set Cover
    beckhams_razor_prots = [r for rg in H.prots() if rg not in HMC for r in rg]
    H.remove_nodes_from([rg for rg in H.prots() if rg not in HMC]) # after that step the drawing will not include small red dots =)
    HF.form_groups(merging_merged=True)# removal of proteins might leave some peptide groups attributed to precisely the same proteins groups
    return H

HE = get_graph(E)
HF = get_graph(F)

# NOW: FINALLY THE BLOODY REPORT
# we also need to merge H nodes again???


pep2pepgroup = {p:pg for pg in H.peps() for p in pg}


# Getting intensities:
D = D.set_index('pep')
D2 = D.loc[pep2pepgroup] #reducing D to peps in R

D.loc["QVIPIIGK"]

# WTF ist das?

D2_pgr = D2.groupby(p2pgr)
I_pgr = D2_pgr[I_cols].sum()# for maximal intensity

p2pgr

# reporting
for rgr in R.prots():
    if len(rgr) > 2:
        break

m = R.subgraph(nx.node_connected_component(R, rgr))
m.draw(with_labels=True)




it = (cc for cc in H.components() if len(cc) > 2)
cc = next(it)
cc.draw(node_size=[40 if n in HMC else 10 for n in cc])







# R = ProtPepGraph((r,p) for pH in H.peps() if H.degree(pH)==1
#                   for r in H[pH] for p in H[r])
# I = ProtPepGraph((r,p) for r,p in H.prot_pep_pairs() if p not in R)
# for e in I.AB():
#     R.add_AB_edge(*e)
# J = ProtPepGraph((r,p) for r,p in H.prot_pep_pairs() if r in I)
# R = R.form_groups(merging_merged=True)
# # return H, R, I, J
# H, R, I, J = simplify(G)

# Facilitate the task by including the supported edges.
# supported = BiGraph((a,b) for c in G.B() if G.degree(c)==1
#                     for a in G[c] for b in G[a])
# supported.draw(with_labels=True)
# unsupported = BiGraph((a,b) for a,b in G.AB() if b not in supported)
# unsupported.draw(with_labels=True)
# input for MSC, if not a cycle



# something is still wrong here: one protein group should have at most one peptide group with deg=1.

I_pgr.loc[R[rgr]].sum() # maximal intensity per run for a protein group
I_pgr.loc[(p for p in R[rgr] if R.degree(p) == 1)].sum() # minimal
# how can we have two peptide groups at this point with deg 1? We don't now.


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

