# From Peptidome to Proteome via Mainzer Dome

Our software aims at infering protein levels from observed peptide levels in a mass spectrometry experiment.

Example use:

```python
import furious_fastas as ff

from pep2prot.graph_ops import get_minimal_protein_group_coverage
from pep2prot.readers import read_ms2rescore_peptide_report

fastas_path = "data/Human_2024_02_16_UniProt_Taxon9606_Reviewed_20434entries_contaminant_tenzer.fasta"
ms2rescore_input = "data/results.sage.ms2rescore.mokapot.peptides.txt"

fastas = [(f.header.split(" ", 1)[0][1:], f.sequence) for f in ff.fastas(fastas_path)]
peptide_report = read_ms2rescore_peptide_report(ms2rescore_input)

covering_protein_groups = get_minimal_protein_group_coverage(
    peptides=list(set(peptide_report.raw_sequence)),
    fastas=fastas,
    min_number_of_peptides=3,
)
```

This should result in a pandas dataframe, like so:

```
      prot_idxs  repr_prot_id  peptide_cnt                    accessions  prot_in_group_cnt                                            pep_ids
group
12         [69]            69            3       [sp|A4UGR9|XIRP2_HUMAN]                  1                                  [454, 1832, 2286]
79        [489]           489            3        [sp|O14939|PLD2_HUMAN]                  1                                  [692, 2571, 2694]
96        [537]           537            4       [sp|O15078|CE290_HUMAN]                  1                            [398, 1583, 1841, 1880]
131       [738]           738            3       [sp|O43295|SRGP3_HUMAN]                  1                                 [1118, 1629, 2230]
132       [743]           743            3       [sp|O43303|CP110_HUMAN]                  1                                   [271, 859, 2877]
...         ...           ...          ...                           ...                ...                                                ...
2517    [20436]         20436           60  [sp|P02769|ALBU_BOVIN_CONTA]                  1  [202, 339, 394, 414, 596, 618, 629, 730, 905, ...
2519    [20566]         20566            3    [sp|P00761|TRYP_PIG_CONTA]                  1                                  [887, 2869, 3111]
2520    [20567]         20567           23  [sp|P00924|ENO1_YEAST_CONTA]                  1  [3, 21, 51, 63, 65, 70, 147, 251, 336, 646, 79...
2521    [20568]         20568           23  [sp|P00925|ENO2_YEAST_CONTA]                  1  [21, 63, 110, 147, 251, 367, 598, 646, 668, 68...
2522    [20569]         20569            9  [sp|P49065|ALBU_RABIT_CONTA]                  1  [618, 1710, 1816, 1893, 2334, 2525, 2573, 2891...
```

# Installation

On Linux.
```
python3 -m venv ve_pep2prot
source ve_pep2prot/bin/activate
pip install -r requirements.txt
```