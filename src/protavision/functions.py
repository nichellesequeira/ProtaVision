from Bio import SeqIO
from Bio import ExPASy
import pandas as pd
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SwissProt
from Bio import pairwise2


def get_protein_sequence(protein_name):
    """
    Retrieve the protein sequence from the ExPASy database using the provided protein name.

    Parameters
    ----------
    protein_name : str
        The name of the protein in the Swiss-Prot database (e.g., "MYG_HUMAN").

    Returns
    -------
    str
        The amino acid sequence of the protein.

    Raises
    ------
    ValueError
        If the provided protein name is not found in the database or if there is an issue with retrieving the sequence.

    Examples
    --------
    >>> get_protein_sequence("MYG_HUMAN")
    'MADQLTEEQIAEFKEAFSLFDKDGDGTITTKE...'
    """
    handle = ExPASy.get_sprot_raw(protein_name)
    record = SeqIO.read(handle, "swiss")
    return str(record.seq)
 

def compare_sequences(sequence1, sequence2):
    """
    Compare two amino acid sequences and return a DataFrame with differences.

    Parameters
    ----------
    sequence1 : str
        The first amino acid sequence.
    sequence2 : str
        The second amino acid sequence.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the index and the differing amino acids of the sequences.

    Raises
    ------
    ValueError
        If the lengths of the sequences are not equal.
    """
    if len(sequence1) != len(sequence2):
        raise ValueError("The lengths of the sequences are not equal.")    
    
    # Create a list to store details of differences
    differences = []
    n = 0
    
    # Go through the two sequences and find the differences
    for i in range(0, min(len(sequence1), len(sequence2))):
        if sequence1[i] != sequence2[i]:
            # Add difference details to the list
            differences.append([i, sequence1[i], sequence2[i]])
            n += 1

    print(f"The total number of amino acids that don't match is: {n}")
    print(f"The total number of amino acids that match is: {len(sequence1) - n}")
    
    # Create a pandas DataFrame from the difference list
    if n == 0:
        return pd.DataFrame()
    else:
        return pd.DataFrame(differences, columns=["Index", "Sequence 1", "Sequence 2"])





import matplotlib.pyplot as plt
import numpy as np

def proportion_amino_acid(sequence1, sequence2):
    """
    Compare two amino acid sequences and generate a bar chart showing the proportion of each amino acid.

    Parameters
    ----------
    sequence1 : str
        The first amino acid sequence.
    sequence2 : str
        The second amino acid sequence.

    Returns
    -------
    None
    """
 
    analyzed_seq1 = ProteinAnalysis(sequence1)
    analyzed_seq2 = ProteinAnalysis(sequence2)


    counts1 = analyzed_seq1.count_amino_acids()
    counts2 = analyzed_seq2.count_amino_acids()

  
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

 
    ax1 = axes[0]
    ax1.bar(counts1.keys(), counts1.values(), color='blue')
    ax1.set_title('Proportions of amino acids - Sequence 1')
    ax1.set_xlabel('Amino Acid')
    ax1.set_ylabel('Number')

    ax2 = axes[1]
    ax2.bar(counts2.keys(), counts2.values(), color='red')
    ax2.set_title('Proportions of amino acids - Sequence 2')
    ax2.set_xlabel('Amino Acid')
    ax2.set_ylabel('Number')

 
    plt.tight_layout()
    plt.show()








def count_conservative_substitutions(sequence1, sequence2):
    """
    Count the number of conservative substitutions between two protein sequences.

    Conservative substitutions are defined based on a predefined conservation matrix,
    which lists amino acids that can be substituted for each other without significantly 
    affecting the protein's function.

    Args:
        sequence1 (str): The first protein sequence.
        sequence2 (str): The second protein sequence.

    Returns:
        int: The number of conservative substitutions between the two sequences.

    Raises:
        ValueError: If the input sequences are not of the same length.

    Example:
        >>> seq1 = "ANQHRQ"
        >>> seq2 = "LMANPR"
        >>> count_conservative_substitutions(seq1, seq2)
        2

    """
    # Define the conservation matrix directly in the functio
    conservation_matrix = {
        'C': {'A'},  # Cysteine
        'S': {'T', 'A', 'G', 'N', 'D', 'E', 'Q', 'K'},  # Serine
        'T': {'S', 'A', 'N', 'V'},  # Threonine
        'P': {},  # Proline
        'A': {'C', 'S', 'T', 'G', 'V'},  # Alanine
        'G': {'S', 'A', 'N', },  # Glycine
        'N': {'S', 'T', 'G', 'D', 'E', 'Q', 'H', 'R', 'K'},  # Asparagine
        'D': {'S', 'N', 'E', 'Q'},  # Aspartic acid
        'E': {'S', 'N', 'D', 'Q', 'H', 'R', 'K'},  # Glutamic Acid
        'Q': {'S', 'N', 'D', 'E', 'H', 'R', 'K', 'M'},  # Glutamine
        'H': {'N', 'E', 'Q', 'R', 'Y'},  # Histidine
        'R': {'N', 'E', 'Q', 'H', 'K'},  # Arginine
        'K': {'S', 'N', 'E', 'Q', 'R'},  # Lysine
        'M': {'Q', 'I', 'L', 'V', 'F'},  # Methionine
        'I': {'M', 'L', 'V', 'F'},  # Isoleucine
        'L': {'M', 'I', 'V', 'F'},  # Leucine
        'V': {'T', 'A', 'M', 'I', 'L'},  # Valine
        'F': {'M', 'I', 'L', 'Y', 'W'},  # Phenylalanine
        'Y': {'H', 'F', 'W'},  # Tyrosine
        'W': {'F', 'Y'},  # Tryptophan
    }

    # Make sure both sequences have the same length
    if len(sequence1) != len(sequence2):
        raise ValueError("The sequences have to be of same length!")

    # Number of conservative substitutions
    conserv_substitutions = 0

    # Run both sequences simultaneously and compare amino acids
    for aa1, aa2 in zip(sequence1, sequence2):
        # Check whether amino acids can be conservatively substituted
        if aa2 in conservation_matrix.get(aa1, set()):
            conserv_substitutions += 1

    return conserv_substitutions



def uniprot_to_pdb(uniprot_id):
    """
    Retrieve the PDB code associated with a given UniProt ID.

    This function searches for UniProt information for the given UniProt ID using the ExPASy database,
    extracts the PDB codes associated with the UniProt sequence, and returns the first PDB code found.

    Args:
        uniprot_id (str): The UniProt ID of the protein.

    Returns:
        str or None: The first PDB code associated with the given UniProt ID, or None if no PDB code is found.

    Example:
        >>> uniprot_to_pdb("MYG_HUMAN")
        '3RGK'

    Raises:
        Exception: If there is an issue retrieving or parsing the UniProt record.

    """

    # Search for UniProt information for the given UniProt ID
    handle = ExPASy.get_sprot_raw(uniprot_id)
    record = SwissProt.read(handle)

    # Retrieve the PDB code associated with the UniProt sequence
    pdb_ids = [ref[1] for ref in record.cross_references if ref[0] == 'PDB']
    if pdb_ids:
        # Select the first PDB code from the list of PDB codes associated with the UniProt protein.
        pdb_id = pdb_ids[0]
        return pdb_id
    else:
        return None



def calculate_alignment_details(sequence1, sequence2):
    """
    Calculate the alignment details between two sequences.
    
    Parameters:
        sequence1 (str): The first sequence.
        sequence2 (str): The second sequence.
    
    Returns:
        tuple: A tuple containing the aligned sequences, start and end positions, and number of matches.
    """
    # Perform global alignment
    alignments = pairwise2.align.globalxx(sequence1, sequence2)
    
    # Extract alignment details from the first alignment (assumes only one alignment is found)
    alignment = alignments[0]
    sequence1_aligned, sequence2_aligned, _, begin, end = alignment
    
    # Count matches
    num_matches = sum(a == b for a, b in zip(sequence1_aligned, sequence2_aligned))
    
    return sequence1_aligned, sequence2_aligned, begin, end, num_matches



def calculate_number_of_gaps(sequence_aligned, sequence_original):
    """
    Calculate the number of gaps in the alignment.
    
    Parameters:
        sequence_aligned (str): The aligned sequence.
        sequence_original (str): The original unaligned sequence.
    
    Returns:
        int: The number of gaps in the alignment.
    """
    num_gaps = sequence_aligned.count("-")
    return num_gaps




def count_matches_with_gap(sequence1, sequence2, gap_length, position):
    """
    Count the number of matches between two sequences with a gap added.
    
    Parameters:
        sequence1 (str): The first sequence.
        sequence2 (str): The second sequence.
        gap_length (int): The length of the gap to be inserted.
        position (int): The position at which to insert the gap.
    
    Returns:
        int: The number of matches between the sequences with the gap.
    """
    # Convert the first sequence to a list to allow insertion
    new_seq = list(sequence1)
    
    # Insert the gap into the new sequence
    for i in range(gap_length):
        new_seq.insert(position + i, "X")
    
    # Determine the length of the shorter sequence
    min_length = min(len(new_seq), len(sequence2))
    
    # Count the number of matches
    count = sum(1 for i in range(min_length) if new_seq[i] == sequence2[i])
    
    return count

def count_amino_acids(sequence):
    """
    Count the number of amino acids in different categories within a given sequence.

    This function categorizes amino acids into hydrophobic, hydrophilic, acidic, and basic groups,
    then counts the occurrences of each category in the provided amino acid sequence.

    Args:
        sequence (str): The amino acid sequence to be analyzed.

    Returns:
        dict: A dictionary with the counts of amino acids in each category:
            - "hydrophobes": Hydrophobic amino acids (A, V, I, L, M, F, Y, W)
            - "hydrophiles": Hydrophilic amino acids (N, C, Q, S, T)
            - "acides": Acidic amino acids (D, E)
            - "bases": Basic amino acids (K, R, H)

    Example:
        >>> count_amino_acids("MVHLTPEEK")
        {'hydrophobics': 3, 'hydrophiles': 1, 'acids': 2, 'bases': 2}
    """
    # Definition of amino acid categories
    hydrophobics = "AVILMFYW"
    hydrophiles = "NCQST"
    acids = "DE"
    bases = "KRH"

    # Counter initialization
    counts = {
        "hydrophobics": 0,
        "hydrophiles": 0,
        "acids": 0,
        "bases": 0
    }

    # Amino acid sequence and count
    for aa in sequence.upper():
        if aa in hydrophobics:
            counts["hydrophobics"] += 1
        elif aa in hydrophiles:
            counts["hydrophiles"] += 1
        elif aa in acids:
            counts["acids"] += 1
        elif aa in bases:
            counts["bases"] += 1

    return counts
