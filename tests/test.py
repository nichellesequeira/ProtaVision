import sys
sys.path.insert(0, '../src/protavision')
import pytest
import functions
from Bio import ExPASy, SeqIO

from functions import get_protein_sequence
def get_protein_sequence_success():
    protein_name = 'MYG_HUMAN'
    expected_sequence = (
        'MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLK'
        'KHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGAD'
        'AQGAMNKALELFRKDMASNYKELGFQG'
    )
    
    result = get_protein_sequence(protein_name)
    assert result == expected_sequence


def get_protein_sequence_no_uniprot_code():
    protein_name = 'MYG_CAT'
    with pytest.raises(ValueError):
        get_protein_sequence(protein_name)
        


def get_protein_sequence_success():
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

def count_conservative_substitutions_none():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "ACDEFGHIKLMNPQRSTVWY"
    assert count_conservative_substitutions(sequence1, sequence2) == 0

def count_conservative_substitutions_diff_len():
    sequence1 = "ACDEFGHIKLMNPQRSTVWY"
    sequence2 = "BCDEFGHIJKLMNOPQRSTUVWXYZ"
    with pytest.raises(ValueError):
        count_conservative_substitutions(sequence1, sequence2)

def count_conservative_substitutions_value2():
    sequence1 = "ANQHRQ"
    sequence2=  "LMANPR"
    assert count_conservative_substitutions(sequence1, sequence2) == 2




    






