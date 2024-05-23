def main() :
    from Bio import SeqIO
    from Bio import ExPASy
    import sys
    import pandas as pd
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    sys.path.append("../src")
    from protavision import functions
    import matplotlib.pyplot as plt
    
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
    
    
    # Print the molecular weight using an f-string
    
    print(f"The molecular weight of your first amino acid sequence is: {molecular_weight1} g/mol")
    print()
    print(f"The molecular weight of your second amino acid sequence is: {molecular_weight2} g/mol")

    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    functions.proportion_amino_acid(sequence1, sequence2)

    num_conserv_substitutions = functions.count_conservative_substitutions(sequence1, sequence2)
    print("The number of conservative substitutions :", num_conserv_substitutions)

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

    import py3Dmol #install package py3Dmol is required
    
    view = py3Dmol.view(query=pdb_id1)
    
    view.setStyle({'cartoon':{'color':'spectrum'}})
    
    view

    import py3Dmol #install package py3Dmol is required
    
    view = py3Dmol.view(query=pdb_id2)
    
    view.setStyle({'cartoon':{'color':'spectrum'}})
    
    view

    from Bio import pairwise2
    
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

    import os
    import ast
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Determine the path to the aminoacids.txt file
    current_dir = os.path.dirname(os.path.abspath('__file__'))
    project_root = os.path.dirname(current_dir)
    file_path = os.path.join(project_root, 'data', 'aminoacids.txt')
    
    # Read the amino acids data from the file
    with open(file_path, 'r') as file:
        data = file.read().strip()
        amino_acids = ast.literal_eval(data)
    
    # Create a DataFrame
    df = pd.DataFrame(amino_acids)
    
    # Display the DataFrame
    print("Amino Acid Data:")
    display(df)
    
    counts_1rsprotein = functions.count_amino_acids(sequence1)
    counts_2ndprotein = functions.count_amino_acids(sequence2)
    print(f"Séquence : {sequence1}")
    print("Number of hydrophobic amino acids :", counts_1rsprotein["hydrophobes"])
    print("Number of hydrophilic amino acids  :", counts_1rsprotein["hydrophiles"])
    print("Number of acidic amino acids :", counts_1rsprotein["acides"])
    print("Number of basic amino acids  :", counts_1rsprotein["basiques"])

    return

