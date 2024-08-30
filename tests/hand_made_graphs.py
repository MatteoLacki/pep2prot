from pep2prot.graph_ops2 import get_minimal_protein_group_coverage


def test_inference_works():
    proteins = [
        (
            ">t|A",
            "BRAAABBBZWE",
        ),
        (
            ">t|B",
            "DAAAABBBADA",
        ),
        (
            ">t|C",
            "WXVBBBASC",
        ),
        (
            ">t|D",
            "ABCBBBCCCDDDXYZ",
        ),
        (">t|E", "ZXCVBNM"),
        (">t|F", "VBNMZXCV"),
    ]
    peptides = ["AAA", "BBB", "CCC", "DDD", "EEE", "FFF"]

    covering_protein_groups = get_minimal_protein_group_coverage(
        peptides=peptides,
        fastas=proteins,
        min_number_of_peptides=1,
    )

    expected_covered_peps_count = len(
        set(
            pep
            for pep in peptides
            if any(pep in prot for prot_header, prot in proteins)
        )
    )
    assert expected_covered_peps_count == len(
        set(
            pep
            for covered_peps in covering_protein_groups.pep_ids
            for pep in covered_peps
        )
    ), "Some peptides were not covered."

    expected_proteins = set(
        prot for prots in covering_protein_groups.accessions for prot in prots
    )

    for prot in (">t|C", ">t|E", ">t|F"):
        assert prot not in expected_proteins
