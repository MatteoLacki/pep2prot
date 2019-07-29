%load_ext autoreload
%autoreload 2

from pathlib import Path
from pep2prot.read import read_fastas
from pep2prot.read import read_isoquant_peptide_report

from furious_fastas.fastas import UniprotFastas


def get_fastas(path):
    u = UniprotFastas()
    u.read(path)
    return u

#typical fasta
path0 = "/home/matteo/Projects/pep2prot/pep2prot/data/mouse.fasta"
fastas0 = get_fastas(path0)

fastas0


data = Path(r"/home/matteo/Projects/pep2prot/tests/sabine3")
path1 = data/'20180515_up_mouse_reviewed_16970entries_and_171contaminants_reverese.fas'
fastas1 = get_fastas(path1) 
path2 = data/'up_Mus_musculus_entries20431_reviewed_reverese.fas'
fastas2 = get_fastas(path2) 



[f for f in fastas0 if 'REVERSE' in f.header]
[f for f in fastas1 if 'REVERSE' in f.header]
[f for f in fastas2 if 'REVERSE' in f.header]

fastas0[0]

[f for f in fastas1 if "Q66JV4" in str(f)]

def parse_fasta_header(f):
    return (' '.join(f.header.split(' ')[2:]), str(f), f.header.split(' ')[1]) 


parse_fasta_header(fastas0[0])
fastas0[0].header

fastas1[0]
fastas2[0]

# open another Sabine's project and try to identify the proteins.
X0 = read_fastas(path0)
X1 = read_fastas(path1)

fastas = UniprotFastas()
fastas.read(path1)
fastas[-100:-90]



import pandas as pd


def parse_fasta_header2(f):
    return (' '.join(f.header.split(' ')[1:]), str(f), f.header.split('|')[2].split(' ')[0]) 

fastas_df = pd.DataFrame.from_records(parse_fasta_header2(f) for f in fastas)
fastas_df.columns = ['description', 'prot_seq', 'prot']
fastas_df = fastas_df.set_index('prot')
fastas_df['seq_len'] = fastas_df.prot_seq.map(len)



D = read_isoquant_peptide_report(data/'pep_rep.csv')
D["prots"] = D.pre_homology_entries.str.split(',').apply(frozenset)

all_prots = {r for rg in D.prots for r in rg}
len(all_prots)

"2A5G_MOUSE" in fastas_df.index
r in fastas_df.index for r in all_prots

strange_prots = [r for r in all_prots if r not in fastas_df.index]
'REVERSE2433'

[r for r in fastas_df.index if 'REVERSE' in r]

