import pytest
from protavision.functions import get_protein_sequence, compare_sequences, proportion_amino_acid, count_conservative_substitutions, uniprot_to_pdb, calculate_alignment_details, calculate_number_of_gaps, count_matches_with_gap, count_amino_acids
import matplotlib.pyplot as plt

def test_get_protein_sequence_success_myoh():
    protein_name = 'MYG_HUMAN'
    expected_sequence = (
        'MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLK'
        'KHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGAD'
        'AQGAMNKALELFRKDMASNYKELGFQG'
    )
    
    result = get_protein_sequence(protein_name)
    assert result == expected_sequence


def test_get_protein_sequence_no_uniprot_code():
    protein_name = 'MYG_CAT'
    with pytest.raises(ValueError):
        get_protein_sequence(protein_name)

def test_get_protein_sequence_success_hemoh():
    protein_name='HBB_HUMAN'
    expected_sequence = (
    'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK'
    'VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG'
    'KEFTPPVQAAYQKVVAGVANALAHKYH'
    )
    result = get_protein_sequence(protein_name)
    assert result == expected_sequence


def test_compare_sequences_same_sequences():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDEFGHIKLMNPQRSTVWY"
    result = compare_sequences(sequence1, sequence2)
    assert len(result) == 0

def test_compare_sequences_different_sequences():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "BCDEFGHIJKLMNOPQRSTU"
    result = compare_sequences(sequence1, sequence2)
    assert len(result) == 13

def test_compare_sequences_different_lengths():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDEFGHIJKLMNOPQRSTUVWXYZ"
    with pytest.raises(ValueError):
        compare_sequences(sequence1, sequence2)

def test_proportion_amino_acid_valid_sequences():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "BCDEFGHIJKLMNOPQRSTUVWXYZ"
    try:
        proportion_amino_acid(sequence1, sequence2)
    except Exception as e:
        pytest.fail(f"Unexpected exception: {e}")

def test_proportion_amino_acid_valid_sequences_with_same_amino_acid():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "AAAAAAAAAAAAAAAAAAAA"
    try:
        proportion_amino_acid(sequence1, sequence2)
    except Exception as e:
        pytest.fail(f"Unexpected exception: {e}")

def test_count_conservative_substitutions_none():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDEFGHIKLMNPQRSTVWY"
    assert count_conservative_substitutions(sequence1, sequence2) == 0

def test_count_conservative_substitutions_diff_len():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "BCDEFGHIJKLMNOPQRSTUVWXYZ"
    with pytest.raises(ValueError):
        count_conservative_substitutions(sequence1, sequence2)

def test_count_conservative_substitutions_value2():
    sequence1 = "ANQHRQ"
    sequence2 = "LMANPR"
    assert count_conservative_substitutions(sequence1, sequence2) == 2


def test_uniprot_to_pdb_valid_id_mygh():
    uniprot_id = "MYG_HUMAN"
    expected_pdb_id = "3RGK"  
    try:
        pdb_id = uniprot_to_pdb(uniprot_id)
        assert pdb_id == expected_pdb_id
    except Exception as e:
        pytest.fail(f"Unexpected exception: {e}")


def test_calculate_alignment_details_perfect_match():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDEFGHIKLMNPQRSTVWY"
    result = calculate_alignment_details(sequence1, sequence2)
    sequence1_aligned, sequence2_aligned, begin, end, num_matches = result
    
    assert sequence1_aligned == sequence1
    assert sequence2_aligned == sequence2
    assert begin == 0
    assert end == len(sequence1)
    assert num_matches == len(sequence1)

def test_calculate_alignment_details_partial_match():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDEXFGHIKLMNPQRSTVWY"
    result = calculate_alignment_details(sequence1, sequence2)
    sequence1_aligned, sequence2_aligned, begin, end, num_matches = result
    
    assert sequence1_aligned == "ACDE-FGHIKLMNPQRSTVWY"
    assert sequence2_aligned == "ACDEXFGHIKLMNPQRSTVWY"
    assert begin == 0
    assert end == len(sequence2)
    assert num_matches == 20


def test_calculate_alignment_details_different_lengths():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"
    result = calculate_alignment_details(sequence1, sequence2)
    sequence1_aligned, sequence2_aligned, begin, end, num_matches = result
    
    expected_sequence1_aligned = sequence1 + "-" * (len(sequence2) - len(sequence1))
    
    assert sequence1_aligned == expected_sequence1_aligned
    assert sequence2_aligned == sequence2
    assert begin == 0
    assert end == len(sequence2)
    assert num_matches == len(sequence1)

def test_calculate_number_of_gaps_success():
    sequence_aligned = "A-CDEFGHIKL-MNPQRSTVWY"
    sequence_original = "ACDEFGHIKLMNPQRSTVWY"
    num_gaps = calculate_number_of_gaps(sequence_aligned, sequence_original)
    expected_num_gaps = 2
    assert num_gaps == expected_num_gaps

def test_calculate_number_of_gaps_no_gaps():
    sequence_aligned = "ACDEFGHIKLMNPQRSTVWY"
    sequence_original = "ACDEFGHIKLMNPQRSTVWY"
    num_gaps = calculate_number_of_gaps(sequence_aligned, sequence_original)
    expected_num_gaps = 0
    assert num_gaps == expected_num_gaps


def test_count_matches_with_gap_no_gap():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDEFGHIKLMNPQRSTVWY"
    gap_length = 0
    position = 0
    expected_matches = len(sequence1)
    assert count_matches_with_gap(sequence1, sequence2, gap_length, position) == expected_matches

def test_count_matches_with_gap_inserted_gap():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDEXFGHIKLMNPQRSTVWY"
    gap_length = 1
    position = 4
    expected_matches = len(sequence2)
    assert count_matches_with_gap(sequence1, sequence2, gap_length, position) == expected_matches

def test_count_matches_with_gap_multiple_gaps():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDXXFGHIKLMNPQRSTVWY"
    gap_length = 1
    position = 3
    expected_matches = len(sequence2) - 1 
    assert count_matches_with_gap(sequence1, sequence2, gap_length, position) == expected_matches


def test_count_matches_with_gap_start_position():
    sequence1 = "CDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDEFGHIKLMNPQRSTVWY"
    gap_length = 1
    position = 0
    expected_matches = len(sequence1)
    assert count_matches_with_gap(sequence1, sequence2, gap_length, position) == expected_matches



def test_count_amino_acids_success():
    sequence = "MVHLTPEEK"
    counts = count_amino_acids(sequence)
    expected_counts = {
        "hydrophobics": 3,
        "hydrophiles": 1,
        "acids": 2,
        "bases":2
    }
    assert counts == expected_counts


