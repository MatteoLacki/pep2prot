%load_ext autoreload
%autoreload 2

import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import networkx as nx
import numpy as np
from pathlib import Path
import difflib
import pandas as pd
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_colwidth', 40)#display whole column without truncation

from furious_fastas.fastas import UniprotFastas
from furious_fastas.contaminants import uniprot_contaminants
from aa2atom import aa2atom, atom2mass
from aa2atom.aa2atom import UnknownAminoAcid

from pep2prot.graphs import ProtPepGraph, BiGraph


path = r"~/Projects/pep2prot/pep2prot/data"

path = Path(path).expanduser()
D = pd.read_csv(path/'peptide_report.csv', encoding = "ISO-8859-1")
fastas = UniprotFastas()
fastas.read(path/'mouse.fasta')

min_pepNo_per_prot = 2

# preprocessing
D.modifier = D.modifier.fillna('')
I_cols = [c for c in D.columns if "intensity in" in c]
D[I_cols] = D[I_cols].fillna(0)
D.rename(columns={'pre_homology_entries':'prots'}, inplace=True)
D['pep'] = np.where(D['modifier'], D['sequence'] + "_" + D['modifier'], D['sequence'])
pep_prots = ['pep','prots']
D.prots = D.prots.str.split(',')
accession2seq = {f.header.split(' ')[1]: str(f) for f in fastas}
obs_accessions = {r for rg in D.prots for r in rg}
assert obs_accessions <= set(accession2seq)


# Fun begins here:   peptide BLUE   protein: RED
G = ProtPepGraph((r,p) for rs, p in zip(D.prots, D.pep) for r in rs)
prots_without_enough_peps = [r for r in G.prots() if G.degree(r) < min_pepNo_per_prot]
G.remove_nodes_from(prots_without_enough_peps)

#TODO: add to ProtPepGraph methods.
def simplify(G):
    H = G.form_groups()
    R = ProtPepGraph((r,p) for pH in H.peps() if H.degree(pH)==1
                           for r in H[pH] for p in H[r]) # result graph
    I = ProtPepGraph((r,p) for r,p in H.prot_pep_pairs() if p not in R) # problematic assignments
    for e in I.AB():
        R.add_AB_edge(*e)
    J = ProtPepGraph((r,p) for r,p in H.prot_pep_pairs() if r in I)
    R = R.form_groups(merging_merged=True)
    return H, R, I, J

# def test_cycles_1():
#     """Test """
#     G = ProtPepGraph([('A',1), ('A',2), ('A',4), ('B',1), ('B',3), ('C',3), ('C',2), ('Z',4), ('Z',5)])
#     return simplify(G)
# H,R,I,J = test_cycles_1()

# def test_cycles_2():
#     """Test """
#     G = ProtPepGraph([('A',1), ('A',2), ('B',1), ('B',3), ('C',3), ('C',2)])
#     return simplify(G)
# H,R,I,J = test_cycles_2()
# H.draw(with_labels=True)
# R.draw(with_labels=True)
# I.draw(with_labels=True)
H, R, I, J = simplify(G)
I.draw(with_labels=True)











# Getting intensities:
D = D.set_index('pep')
D2 = D.loc[(r for rg in R.peps() for r in rg)] #reducing D to peps in R
p2pgr = {p:pg for pg in R.peps() for p in pg}
D2_pgr = D2.groupby(p2pgr)
I_pgr = D2_pgr[I_cols].sum()

# reporting
for rgr in R.prots():
    if len(rgr) > 2:
        break

m = R.subgraph(nx.node_connected_component(R, rgr))
m.draw(with_labels=True)
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

