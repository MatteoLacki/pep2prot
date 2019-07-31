%load_ext autoreload
%autoreload 2

from pathlib import Path
import pandas as pd

from furious_fastas.fastas import Fastas

from pep2prot.read import read_fastas
from pep2prot.read import read_isoquant_peptide_report


def get_fastas(path):
    u = Fastas()
    u.read(path)
    return u


#typical fasta
data_path = Path("/home/matteo/Projects/pep2prot/pep2prot/data")
mouse = Fastas()
path0 = data_path/'mouse.fasta'
mouse.read(path0)

data = Path(r"/home/matteo/Projects/pep2prot/tests/sabine3")
path1 = data/'20180515_up_mouse_reviewed_16970entries_and_171contaminants_reverese.fas'
fastas1 = get_fastas(path1) 
path2 = data/'up_Mus_musculus_entries20431_reviewed_reverese.fas'
fastas2 = get_fastas(path2) 
fastas2.fasta_types()
fastas2[-1]



path = path1
fastas = Fastas()
fastas.read(path)
fastas.same_fasta_types()


fastas_df = pd.DataFrame.from_records((f.description, f.sequence, f.entry) for f in fastas)
fastas_df.columns = ['description', 'prot_seq', 'prot']
fastas_df = fastas_df.set_index('prot')
fastas_df['seq_len'] = fastas_df.prot_seq.map(len)

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

