# Python Portfolio: Task 2

# The code reads the DNA sequences from the FASTA file into a dictionary 
def read_fasta(input_file):
    sequences = {}         
    with open(input_file, 'r') as file:   # this code opens the file 
        sequence_name = ''  
        sequence_data = ''  # this stores the current DNA sequence

        # I created a loop over each line in the FASTA file
        for line in file:
            line = line.strip()    # removes the newline character to avoid /n in the sequence

            # If the line starts with '>', this is a new sequence header
            if line.startswith('>'):
                if sequence_name != '':
                    sequences[sequence_name] = sequence_data

                # Codes to start a new sequence and removes '>' and extra spaces
                sequence_name = line[1:].strip()
                sequence_data = ''      # resets the DNA string

            else:
                sequence_data += line.upper()

        # After the loop finishes, the last sequence is saved
        if sequence_name != '':
            sequences[sequence_name] = sequence_data

    return sequences


# Read the sequences from the FASTA file
# (Changes the filename if your file is called something else)
gene_sequences = read_fasta('Python_task_2_sequences.fa')


# 2. Define the codon tables. Assignment requires three different genetic codes as not all sequences use the standard code.
# Therefore, 3 different codon tables should be defined. 
# Standard code for sequence 1, alterative yeast nuclear code for sequence 2 and vertebrate mitochondrial code for sequence 3.

# Standard genetic code (provided in assignment) for sequence 1
standard_codon_table = {
    # Phenylalanine and Leucine
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    # Leucine
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    # Isoleucine and Methionine (start)
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    # Valine
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    # Serine
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    # Proline
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    # Threonine
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    # Alanine
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    # Tyrosine and Stop
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    # Histidine and Glutamine
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    # Asparagine and Lysine
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    # Aspartic acid and Glutamic acid
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    # Cysteine, Stop, Tryptophan
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    # Arginine
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    # Serine and Arginine
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    # Glycine
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# The alternative yeast nuclear code for sequence 2
# From the NCBI link provided  the main difference was that the CTG codes for Serine (S) instead of Leucine (L)
alt_yeast_nuclear_table = standard_codon_table.copy()
alt_yeast_nuclear_table['CTG'] = 'S'

# The vertebrate mitochondrial code for sequence 3
# From the NCBI link, the differences from the standard code are ATA -> M (Met) instead of I (Ile), TGA -> W (Trp) instead of stop and AGA, AGG -> stop (*) instead of R (Arg)
vertebrate_mito_table = standard_codon_table.copy()
vertebrate_mito_table['ATA'] = 'M'
vertebrate_mito_table['TGA'] = 'W'
vertebrate_mito_table['AGA'] = '*'
vertebrate_mito_table['AGG'] = '*'

# A dictionary was created to match each sequence to the correct codon table
sequence_to_codon_table = {
    'Sequence 1': standard_codon_table,
    'Sequence 2': alt_yeast_nuclear_table,
    'Sequence 3': vertebrate_mito_table
}


# The amino acid property groups (from Table 2 in the assignment)
Hydrophobic = set(['G', 'A', 'V', 'L', 'I', 'M', 'P', 'F', 'W'])
Hydrophilic = set(['S', 'T', 'Y', 'N', 'Q', 'C', 'D', 'E', 'K', 'R', 'H'])
Negatively_Charged = set(['D', 'E'])
Positively_Charged = set(['K', 'R', 'H'])
Uncharged = set(['G', 'A', 'V', 'L', 'I', 'M', 'P', 'F', 'W','S', 'T', 'Y', 'N', 'Q', 'C'])


# Translates DNA into protein using the codon table
# This code converts the DNA sequence into an amino acid chain by reading in groups of three bases (codons),matches to the amino acids and stops at the first stop codon
def translate_gene(dna_sequence, codon_table):
    start_index = dna_sequence.find('ATG')
# If no ATG is found, the code returns an empty string
    if start_index == -1:
        return ''
    amino_acids = '' 
# A loop goes over the sequence in steps of 3 bases
    for position in range(start_index, len(dna_sequence), 3):
        codon = dna_sequence[position:position + 3]
# If fewer than 3 bases are left, stop the loop
        if len(codon) < 3:
            break
# Convert the codon to amino acid using the codon table
        amino_acid = codon_table.get(codon, 'X')  # The letter X is used as an unknown codon
# The loop stops if a stop codon is reached
        if amino_acid == '*':
            break
# Add the amino acid to the protein chain list 
        amino_acids += amino_acid

    return amino_acids


# The amino acid property is counted and then returned in a dictionary
def count_amino_acid_properties(amino_chain):
    Hydrophobic_Count = 0
    Hydrophilic_Count = 0
    Uncharged_Count = 0
    Positively_Charged_Count = 0
    Negatively_Charged_Count = 0

    # A Loop is used over each amino acid in the protein sequence to count the properties
    for aa in amino_chain:
        if aa in Hydrophobic:
            Hydrophobic_Count += 1
        if aa in Hydrophilic:
            Hydrophilic_Count += 1
        if aa in Uncharged:
            Uncharged_Count += 1
        if aa in Positively_Charged:
            Positively_Charged_Count += 1
        if aa in Negatively_Charged:
            Negatively_Charged_Count += 1

    # Afterwards it is returned to the counts as a dictionary
    return {
        'Hydrophobic': Hydrophobic_Count,
        'Hydrophilic': Hydrophilic_Count,
        'Uncharged': Uncharged_Count,
        'Positively Charged': Positively_Charged_Count,
        'Negatively Charged': Negatively_Charged_Count
    }


# The amino acid chain and characteristics are printed 
for sequence_id, dna_seq in gene_sequences.items():

    print(f'\n--- {sequence_id} ---')

    # Choose the correct codon table for this sequence.
    # If the name is not found,the code returns to the default standard code.
    codon_table = sequence_to_codon_table.get(sequence_id, standard_codon_table)
    protein = translate_gene(dna_seq, codon_table)

    # If no start codon was found, there is no ORF starting at ATG
    if protein == '':
        print('No start codon found. The translation is not possible.')
        continue

    # Prints the protein length, amino acids, the properties as a string
    print(f'\nProtein length: {len(protein)} amino acids\n')
    print('Amino acid chain (1-letter code):')
    print(protein)
    properties = count_amino_acid_properties(protein)
    print('\nAmino acid properties:')
    for property_type, count in properties.items():
        print(f'{property_type}: {count}')
