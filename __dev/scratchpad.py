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


# missing which peptides are covered.


# comment: of course, merging peptides makes absolutely no sense in terms of a set cover. So we are right not to do so.
# what now?
# * Missing some form of reporting: for each peptide get a
# * More adapters for different peptide reports
# * Add the q-value calculations as preprocessing:
#   * we need ion level ones
# *
