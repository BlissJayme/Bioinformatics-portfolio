# Task 4: Genomics analyses; variant density or heterozygosity

# 1. Use the scikit-allel package to read the VCF file.
# Imported allel and matplotlib packages to analyse and plot data from a VCF file
# Set working directory to where the VCF file is located on a personal computer

import matplotlib.pyplot as plt
import allel

# The VCF file is read using the allel.read_vcf() function
eels_file = "eels_chromosomes_1-6.vcf.gz"
eel_chromosomes_vcf_file = allel.read_vcf("eels_chromosomes_1-6.vcf.gz")

# The vcf file results are stored in a dictionary format, allowing access to specific data categories using keys
print(eel_chromosomes_vcf_file.keys())



# 2.  Choose one of the following plots to make:
# Number of variants per 100,000 base pair window across one chromosome



# 3. Create your chosen plot for your chosen chromosome or sample group, using matplotlib and allel
# From the keys, the "variants/CHROM" and "variants/POS" keys are used to extract chromosome names and variant positions
# Variants/CHROM is a variant annotation, this is a chromosome identifer for each variant
# Variants/POS is the position of each variant on the chromosom, using 1 based genomic coordinates
chrom = eel_chromosomes_vcf_file["variants/CHROM"]
pos = eel_chromosomes_vcf_file["variants/POS"]

# Extracting variant positions for a single chromosome, Chr_01, and counting variants in 100,000 base pairs
# Creates a list to store positions for the chosen chromosome
# Loops through all variant positions and checks if they belong to the chosen chromosome
# If true, appends the position to the list and is printed
chrom_name = "Chr_01"  
window_size = 100000  
chosen_positions = []  

for i in range(len(pos)):
    if chrom[i] == chrom_name:
        chosen_positions.append(pos[i])

# The window_count takes the position and a window size as inputs, then returns counts and locations for each window
counts, locations = allel.windowed_count(chosen_positions, size=window_size)

# Similar to above, loops through to produce a print out of the number of variants found in each window
# Extracts the number of variants, start and stop positions from results and prints them
# String format is used to present variants found 
for i in range(len(counts)):
    
    variants = counts[i]
    start_position = locations[i][0]
    stop_position = locations[i][1]

    print(f"{variants} variants found between {start_position} and {stop_position}")

# For each window, the midpoint is calculated and stored in a list
midpoints = []
for loc in locations:
    start = loc[0]
    stop = loc[1]
    midpoint = (start + stop) / 2
    midpoints.append(midpoint)

# A graph is created to display the number of variants per 100,000 base pair window across the chromosome
plt.title("The number of variants per 100,000 base pair window across one chromosome")
plt.plot(midpoints, counts)
plt.xlabel("The position on chromosomes")
plt.ylabel("Number of variants")

plt.show()

# A function is created to calculate variant positions for each chromosome
# An empty is formed to store variant positions for the specific chromosome
# A loop goes through all variant positions, checking if they belong to the chosen chromosome
# If true, the position is added to the list
def variant_pos_chromosomes(pos, chrom, chrom_name, window_size):
    chrom_positions = []  
    for i in range(len(pos)):
        if chrom[i] == chrom_name:  
            chrom_positions.append(pos[i])  

# The windowed_count function is used to count variants in specified window sizes
# An empty list is created to store midpoints of each window
# The string prints the variant counts for the current chromosome
    counts, locations = allel.windowed_count(chrom_positions, size=window_size)

    midpoints = []

    print(f"\nVariant counts for {chrom_name}:")

# Loops through each window to calculate midpoints and print variant counts
# The midpoint is calculated by averaging the start and stop positions of each window
# The midpoint is appended to the midpoints list
# A string is used to print the number of variants found in each window
# Returns the midpoints and counts
    for i in range(len(counts)):
        start = locations[i][0]  
        stop = locations[i][1]
        midpoint = (start + stop) / 2  
        midpoints.append(midpoint)

        print(f"{counts[i]} variants found between {start} and {stop}")  

    return midpoints, counts  



# 4. Convert part of your code into a function, using arguments and return, to calculate the same data for every chromosome or sample group
# List of chromosome names to analyse, count variants in 100,000 base pair windows for each chromosome
chromosome_names = ['Chr_01', 'Chr_02', 'Chr_03', 'Chr_04', 'Chr_05', 'Chr_06']
window_size = 100000

# Loops through each chromosome and calculates variant counts 
for chrom_name in chromosome_names:
    midpoints, counts = variant_pos_chromosomes(pos, chrom, chrom_name, window_size)  



# 5. Create a figure or figures which clearly displays this data
# The figure displays the number of variants per 100,000 base pair window for a selected chromosome
# The midpoints and counts are obtained from the function for the chosen chromosome
plt.figure(figsize=(10, 5))  
plt.tight_layout() 
plt.title(f"Number of variants per 100,000 bp window on {chrom_name}")  
plt.plot(midpoints, counts,)  
plt.xlabel("Position on the chromosomes (bp)") 
plt.ylabel("The number of variants") 

plt.show() 

# Acknowledgement of GAI Contribution: Used to assist with the working directory, tigher layout for the matplots to avoid overlapping and a scan of the code for final errors.  


