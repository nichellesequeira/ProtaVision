def main() :
    from Bio import SeqIO
    from Bio import ExPASy
    import sys
    import pandas as pd
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    sys.path.append("../src")
    from protavision import functions
    import matplotlib.pyplot as plt
    import py3Dmol
    from Bio import pairwise2
    
    protein_name1 = input("Please enter the name of the first protein (example MYG_HUMAN): ")
    
    # Call the function using the module name
    sequence1 = functions.get_protein_sequence(protein_name1)
    print(f"This is the sequence of your first protein: {sequence1}")

    protein_name2 = input("Please enter the name of the second protein (example MYG_MOUSE): ")
    sequence2 = functions.get_protein_sequence(protein_name2)
    print(f"This is the sequence of your second protein; {sequence2}")

    
    df = functions.compare_sequences(sequence1, sequence2)
    df
    
    analyzed_seq1 = ProteinAnalysis(str(sequence1))
    analyzed_seq2 = ProteinAnalysis(str(sequence2))
    
    molecular_weight1 = round(analyzed_seq1.molecular_weight(), 3)
    molecular_weight2 = round(analyzed_seq2.molecular_weight(), 3)
    
    
    
    print(f"The molecular weight of your first amino acid sequence is: {molecular_weight1} g/mol")
    print()
    print(f"The molecular weight of your second amino acid sequence is: {molecular_weight2} g/mol")

    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    functions.proportion_amino_acid(sequence1, sequence2)

    num_conserv_substitutions = functions.count_conservative_substitutions(sequence1, sequence2)
    print("The number of conservative substitutions :", num_conserv_substitutions)

    counts_1rsprotein = functions.count_amino_acids(sequence1)
    print(f"Séquence : {sequence1}")
    print("Number of hydrophobic amino acids :", counts_1rsprotein["hydrophobics"])
    print("Number of hydrophilic amino acids  :", counts_1rsprotein["hydrophiles"])
    print("Number of acidic amino acids :", counts_1rsprotein["acids"])
    print("Number of basic amino acids  :", counts_1rsprotein["bases"])

    counts_2ndprotein = functions.count_amino_acids(sequence2)
    print(f"Séquence : {sequence2}")
    print("Number of hydrophobic amino acids :", counts_2ndprotein["hydrophobics"])
    print("Number of hydrophilic amino acids  :", counts_2ndprotein["hydrophiles"])
    print("Number of acidic amino acids :", counts_2ndprotein["acids"])
    print("Number of basic amino acids  :", counts_2ndprotein["bases"])

    pdb_id1 = functions.uniprot_to_pdb(protein_name1)
    if pdb_id1:
        print(f"The PDB Code associated to the first protein {protein_name1} is : {pdb_id1}")
    else:
        print(f"No PDB Code is associated to this first protein {protein_name1}")
    
    pdb_id2 = functions.uniprot_to_pdb(protein_name2)
    if pdb_id2:
        print(f"The PDB code of the second protein {protein_name2} is : {pdb_id2}")
    else:
        print(f"No PDB code is associated to the second protein {protein_name2}")
 
    
    view = py3Dmol.view(query=pdb_id1)
    
    view.setStyle({'cartoon':{'color':'spectrum'}})
    
    view

    
    view = py3Dmol.view(query=pdb_id2)
    
    view.setStyle({'cartoon':{'color':'spectrum'}})
    
    view

    
    sequence1_aligned, sequence2_aligned, begin, end, num_matches = functions.calculate_alignment_details(sequence1, sequence2)
    num_gaps_sequence1 = functions.calculate_number_of_gaps(sequence1_aligned, sequence1)
    num_gaps_sequence2 = functions.calculate_number_of_gaps(sequence2_aligned, sequence2)
    score = 10*num_matches - 25*num_gaps_sequence1
    
    print("Sequence 1 aligned:", sequence1_aligned)
    print("Sequence 2 aligned:", sequence2_aligned)
    print("Number of matches:", num_matches)
    print("Score:", score)
    print("Number of gaps in Sequence 1:", num_gaps_sequence1)
    print("Number of gaps in Sequence 2:", num_gaps_sequence2)

    gap_length = int(input("Enter the length of the gap: "))
    position = int(input("Enter the position to apply the gap: "))
    
    num_matches = functions.count_matches_with_gap(sequence1, sequence2, gap_length, position)
    print("Number of matches found with the gap method is:", num_matches)

    

    return

