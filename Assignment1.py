'''
Team 8

A program to read a RNA sequence from a fasta file, pre-process it
by removing the start codon, translate it according to the standard
genetic code, and print the resulting amino acid sequence to the screen.

'''
import GenetiCode as gene


# Converts RNA sequences to amino acid sequences and returns the result as a string
def rna2aa(rna):
    # replace the following line with the real code of this function
    # get the amino acid dictionary.
    # This is the way to access dictionary elements.
    # print(gene.code["UUG"])
    # print("Debugger for dict")
    # print(gene.code[rna[3:6]])
    # print("Debugger from rna2aa method")
    # print(rna[0:10])
    not_finished = True
    amino_acid = ''
    index = 0
    while not_finished:
        base_pair = rna[index: index + 3]
        # Check if the length of base pair is 3 or not.
        if len(base_pair) < 3:
            # We have finished reading the sequence.
            not_finished = False
        else:
            # We have to check through the dictionary for matching
            # the rna to amino acid sequence.
            amino_acid = amino_acid + gene.code[base_pair]
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


# Calls methods to read gene file and convert RNA to Amino Acids and prints the results
def main():
    genes = read_file("Assignment1Sequences.txt")

    counter = 0
    for rna in genes:
        print("Sequence " + str(counter) + ": " + rna)
        aa_sequence = rna2aa(rna)
        print("The resulting amino acid sequence for gene " + str(counter) + ": " + aa_sequence + "\n")
        counter += 1


main()
