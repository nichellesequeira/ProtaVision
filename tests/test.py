import pytest
import functions
from functions import get_protein_sequence
from Bio import ExPASy, SeqIO
from unittest.mock import patch, MagicMock

def test1_get_protein_sequence():
    protein_name = 'MYG_HUMAN'
    expected_sequence = (
        'MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLK'
        'KHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGAD'
        'AQGAMNKALELFRKDMASNYKELGFQG'
    )
    
    result = get_protein_sequence(protein_name)
    assert result == expected_sequence


def test2_get_protein_sequence():
    protein_name = 'MYG_CAT'
    with pytest.raises(ValueError):
        get_protein_sequence(protein_name)
        


def test3_get_protein_sequence():
    protein_name='HBB_HUMAN'
    expected_sequence = (
    'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK'
    'VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG'
    'KEFTPPVQAAYQKVVAGVANALAHKYH'
    )
    result = get_protein_sequence(protein_name)
    assert result == expected_sequence

import pytest
from functions import count_conservative_substitutions

def test1_count_conservative_substitutions():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDEFGHIKLMNPQRSTVWY"
    assert count_conservative_substitutions(sequence1, sequence2) == 0

def test2_count_conservative_substitutions():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "BCDEFGHIJKLMNOPQRSTUVWXYZ"
    with pytest.raises(ValueError):
        count_conservative_substitutions(sequence1, sequence2)

def test3_count_conservative_substitutions():
    sequence1 = "ANQHRQ"
    sequence2=  "LMANPR"
    assert count_conservative_substitutions(sequence1, sequence2) == 2




    






