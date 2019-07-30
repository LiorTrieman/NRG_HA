# Open and read a fasta file
# analys with regulatory expressions
# output the following:
"""
a. Marker name - e.g. USDA_SNP0012
b. Marker ID - e.g. i00001Gh
c. Marker allele 1 - e.g. C
d. Marker allele 2 - e.g. T
e. Full marker sequence with allele 1 - e.g. TGCAGAACACAGA...C...AAGTAAAA
f. Full marker sequence with allele 2 - e.g. TGCAGAACACAGA...T...AAGTAAAA
g. Marker length - e.g. 63
h. Position of SNP in the marker (0-based index) - e.g. 54
"""
import re #in order to use regular expressions


f = open("TAMU_SNP63K_69997.fasta", "r")
print(f.read())
