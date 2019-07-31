%load_ext autoreload
%autoreload 2

import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 30)#display whole column without truncation

from pathlib import Path
from collections import Counter

from pep2prot.read import read_isoquant_peptide_report, read_fastas
from pep2prot.graphs import ProtPepGraph, BiGraph, get_peptide_protein_graph
from pep2prot.preprocessing import preprocess_isoquant_peptide_report, get_protein_coverages, complex_cluster_buster, simple_cluster_buster 
from pep2prot.intensities import get_prot_intensities
from pep2prot.postprocessing import summarize_prots, get_stats, prettify_protein_informations, get_full_report


test_data = Path(r"/home/matteo/Projects/pep2prot/tests/hye")
pep_rep_path = test_data/'2019-019 HYE neu_user designed 20190621-083010_peptide_quantification_report.csv'
fastas_path = test_data/'HYE.fasta'
cluster_buster=complex_cluster_buster

D = read_isoquant_peptide_report(pep_rep_path)
D, I_cols = preprocess_isoquant_peptide_report(D)
fastas = read_fastas(fastas_path)
observed_prots = {r for rg in D.prots for r in rg}
assert all(r in fastas.index for r in observed_prots), "It seems that you are using a different set of fastas than the peptide annotation software before. Repent please."
prots     = get_protein_coverages(D, fastas)
uni_cols  = ['peptide_overall_max_score','peptide_fdr_level',
                  'peptide_overall_replication_rate','prots',
                  'pre_homology_accessions','pi','mw']

DD = cluster_buster(D, I_cols, uni_cols) # agg same peptides in various clusters
assert np.all(DD.groupby(['pep','prots']).size() == 1), "Some peptides are still across different clusters."

H, prots_no_peps, peps_no_prots, beckham_prots = get_peptide_protein_graph(DD) 
pep2pepgr      = {p:pg for pg in H.peps() for p in pg}
DDinH          = DD.loc[pep2pepgr] # peps in H: no simple prot-pep pairs, no unnecessary prots?
DDinH['pepgr'] = DDinH.index.map(pep2pepgr)
peps_I         = DDinH[I_cols].groupby(pep2pepgr).sum()# peptide groups intensities
peps_I.isna().any().any()

buggy_weights = weights[weights.isnull().any(axis=1)]
buggy_pep_groups = set(buggy_weights.index.get_level_values('pep'))
buggy_peps = {p for pg in buggy_pep_groups for p in pg}
D.loc[buggy_peps]
buggy_prot_groups = {rg for pg in buggy_pep_groups for rg in H[pg]}
H.subgraph(buggy_pep_groups|buggy_prot_groups).draw(with_labels=True) # looks normal

def gua(rg):
    yield rg
    for p in H[rg]:
        yield p
        for r in H[p]:
            yield r
            for p2 in H[r]:
                yield p2


from pep2prot.graphs import ProtPepGraph

test = ProtPepGraph([('A',0), ('A',1), ('B',1), ('B',2), ('B',4), ('C',2), ('C',3), ('C',4), ('C',5), ('D',4), ('D',5), ('D',6)])

def create_G(pep2prots, min_pepNo_per_prot=2):
    G = ProtPepGraph((r,p) for rs, p in zip(pep2prots.prots, pep2prots.index) for r in rs)
    # removing pairs r-p, where both r and p have no other neighbors
    prots_no_peps = {r for r in G.prots() if G.degree(r) < min_pepNo_per_prot}
    peps_no_prots = {p for r in prots_no_peps for p in G[r] if G.degree(p) == 1} 
    G.remove_nodes_from(prots_no_peps)
    G.remove_nodes_from(peps_no_prots)
    return G

def minimize_graph(G):
    H = G.form_groups()
    HMC = H.greedy_minimal_cover() # Her Majesty's Minimal Set Cover.
    beckham_prots = {r for rg in H.prots() if rg not in HMC for r in rg}
    H.remove_nodes_from([rg for rg in H.prots() if rg not in HMC]) # after that step the drawing will not include small red dots =)
    H = H.form_groups(merging_merged=True)# removal of proteins might leave some peptide groups 
    return H, HMC, beckham_prots

H_test, HMC_test, beckham_prots_test = minimize_graph(test)
H_test.draw(with_labels=True)

H.subgraph(nx.node_connected_component(H, rg)).draw(with_labels=True)
H.subgraph(set(gua(rg))).draw(with_labels=True)


buggy_pepprot_groups = {(pg,rg) for pg in buggy_pep_groups for rg in H[pg]}


pg,rg = next(iter(buggy_pepprot_groups))
weights.xs(pep=pg, prot=rg)

weights.loc[pg]
weights.index
pg in weights.index.get_level_values('pep')
rg in weights.index.get_level_values('prot')
(pg, rg) in weights.index

prots_min_I, prots_I, prots_max_I = get_prot_intensities(H, peps_I)

weights.loc[buggy_pep_groups]
x = weights.xs(pg, level='pep')
x = x.isna().any(axis=1)
rg = next(iter(x[x].index))
D.loc['LAADDFR', I_cols]
# there might be a need for left or right merge somewhere

weights.xs((pg,rg), level=('pep','prot'))
weights
prots_curr_I.loc[rg]
prots_curr_I.loc[]
prots_curr_I.loc[[frozenset({'K1C15_HUMAN', 'K1C15_HUMAN_CONTA'})]]
prots_curr_I.loc[[rg]]
# OK, so there is a bug! some protein group is simply not there. how come???
set(otherpeps2mixprots.get_level_values('prot')) - set(prots_curr_I.index)
H[rg]




set(peps2prots.get_level_values('pep')) - set(otherpeps_I.index)



prot_info      = summarize_prots(H, fastas, prots.pep_coverage)
prots_I_nice   = prettify_protein_informations(prots_I, prot_info)
prots_stats    = get_stats(prots_min_I, prots_I, prots_max_I)

if verbose:
    print('Preparing reports.')
all_prots      = get_full_report(prots_min_I, prots_I, prots_max_I)
all_prots_nice = prettify_protein_informations(all_prots, prot_info)




