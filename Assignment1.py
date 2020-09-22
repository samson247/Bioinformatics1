'''
Team 8

A program to read a RNA sequence from a fasta file, pre-process it
by removing the start codon, translate it according to the standard
genetic code, and print the resulting amino acid sequence to the screen.

'''
import GenetiCode as genetic_code
import KyteDoolittleValues as hydropathy_indices
import HelixPropensities as helix_propensities


# Converts RNA sequences to amino acid sequences and returns the result as a string
def rna2aa(rna):
    not_finished = True
    amino_acid = ''
    index = 0
    while not_finished:
        base_pair = rna[index: index + 3]
        # Check if the length of base pair is 3 or not.
        if len(base_pair) < 3 or base_pair in ["UAA", "UAG", "UGA"]:
            # We have finished reading the sequence.
            not_finished = False
        else:
            # We have to check through the dictionary for matching
            # the rna to amino acid sequence.
            amino_acid = amino_acid + genetic_code.code[base_pair]
        index += 3
    return amino_acid


# Reads the FASTA file and returns the RNA sequences in a list
def read_file(filename):
    with open(filename) as geneFile:
        genes = geneFile.read().splitlines()
    for x in genes:
        if ">Gene" in x:
            genes.remove(x)
    genes = filter(None, genes)
    return genes


# Prints the results of the hydrophobicity scan, including first and last indices, amino acid sequence,
# and average scores for a window when they are greater than 1.6
def print_hydrophobicity_results(results):
    print "Hydrophobicity Results (Score > 1.6 signifies potential transmembrane domain)"
    if not results:
        print "No scores > 1.6 were found"
    else:
        print "| Beg | End |     AA Sequence     | Score |"
        print "-------------------------------------------"
        for result in results:
            output = "| "
            for x in result:
                output += x
                output += " | "
            print output
    print "\n"


# Prints the results of the helix propensity scan, including first and last indices, amino acid sequence,
# and average scores for a window when they are less than 0.4
def print_helix_propensity_results(results):
    print "Helix Propensity Results (Score < 0.4)"
    if not results:
        print "No scores < 0.4 were found"
    else:
        print "| Beg | End |     AA Sequence     | Score |"
        print "-------------------------------------------"
        for result in results:
            output = "| "
            for x in result:
                output += x
                output += " | "
            print output
    print "\n"


# Scans the amino acid sequence for hydrophobicity and helix propensity
def scan_sequence(aa_sequence):
    window_size = 19
    hydropathy_score = 0
    helix_propensity_score = 0
    hydropathy_results = []
    helix_propensity_results = []

    start_index = 0
    while start_index < (len(aa_sequence) - window_size):
        aa_subsequence = aa_sequence[start_index: start_index + window_size]
        for aa in aa_subsequence:
            hydropathy_score += hydropathy_indices.indices[aa]
            helix_propensity_score += helix_propensities.helix_propensities[aa]
        if (hydropathy_score / window_size) > 1.6:
            hydropathy_results.append([str(start_index), str(start_index + (window_size - 1)), aa_subsequence,
                                             str(round(hydropathy_score / window_size, 3))])
        if (helix_propensity_score / window_size) < 0.4:
            helix_propensity_results.append([str(start_index), str(start_index + (window_size - 1)), aa_subsequence,
                                             str(round(helix_propensity_score / window_size, 3))])
        start_index += 1
        hydropathy_score = 0
        helix_propensity_score = 0

    print_hydrophobicity_results(hydropathy_results)
    print_helix_propensity_results(helix_propensity_results)


# Calls methods to read gene file and convert RNA to Amino Acids and prints the results
def main():
    genes = read_file("Assignment1Sequences.txt")

    counter = 1
    for rna in genes:
        #print("Sequence " + str(counter) + ": " + rna)
        aa_sequence = rna2aa(rna)
        print("The resulting amino acid sequence for gene " + str(counter) + ": " + aa_sequence + "\n")
        scan_sequence(aa_sequence)
        counter += 1


main()
