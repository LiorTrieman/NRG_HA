# Open and read a fasta file
# analyse with regular expressions
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
import re  # in order to use regular expressions

with open('TAMU_SNP63K_69997.fasta', 'r') as file: # read the fasta file
    data = file.read().replace('\n', '')
print(type(data))  # view data type
# print(data)   # view the data (should be comment to shorten run-time)

marker_name = re.findall(r'USDA_SNP\d{4}', data)  # a.Marker name
print(marker_name)
marker_ID = re.findall(r'[i]\d{5}\w{2}', data)  # b.Marker ID
marker_alleles = re.findall(r's=\w{1}\S\w{1}', data)  # both allels

# need to separate each allele from the above string
marker_len = len(marker_alleles)

print("len= ", marker_len)
print("typ: ", type(marker_alleles))

# to find  Full marker sequence we should search for :spaceX2, {A,G,C,T, and Y,M,K,R,W,S} and stop at "enter" or "<"

full_marker_seq_raw = re.findall(r'\s\s[A, G, C, T, Y, M, K, R, W, S]*', data)  # Full marker sequence with allele 1
print(full_marker_seq_raw[0:10])
full_marker_seq_raw_len = len(full_marker_seq_raw)
print("len: ", full_marker_seq_raw_len)


def put_allel(text, allele):  # replace one of [Y, M, K, R, W, S] with the alleles
    for ch in ['Y', 'M', 'K', 'R', 'W', 'S']:
        if ch in text:
            text = text.replace(ch, allele) # replace with one of the alleles
    return text


for index in range(0, 2):  # marker_len
    Allele = marker_alleles[index]
    Allele_1 = Allele[2:3]   # c. Marker allele 1
    Allele_2 = Allele[4:5]   # c. Marker allele 2
    full_marker_seq_single_marker = full_marker_seq_raw[index]
    full_marker_seq = full_marker_seq_single_marker[2:]
    # full_marker_seq_allel_1 = full_marker_seq.replace("Y", Allele_1)
    # full_marker_seq_allel_2 = full_marker_seq.replace("Y", Allele_2)
    full_marker_seq_allel_1 = put_allel(full_marker_seq, Allele_1)
    full_marker_seq_allel_2 = put_allel(full_marker_seq, Allele_2)
    print(full_marker_seq_allel_1)
    print(full_marker_seq_allel_2)

'''
full_marker_seq_allele_2 = re.findall(r'[i]\d{5}\w{2}', data)  # Full marker sequence with allele 2
marker_len = re.findall(r'[i]\d{5}\w{2}', data)  # Marker length
snp_position = re.findall(r'[i]\d{5}\w{2}', data)  # Position of SNP in the marker (0-based index)
# 





import bio
help(bio)

#for record in bio.SeqIO.parse("TAMU_SNP63K_69997.fasta", "fasta"):
 #   print(record.id)
'''