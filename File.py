with open('Assignment1Sequences.txt') as geneFile:
    genes = geneFile.read().splitlines()
for x in genes:
    if ">Gene" in x:
        genes.remove(x)
genes = filter(None, genes)
print(genes)
