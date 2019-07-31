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
import csv # to create a new csv file with the new extracted data
with open('TAMU_SNP63K_69997.fasta', 'r') as file: # read the fasta file
    data = file.read().replace('\n', '')
print(type(data))  # view data type
# print(data)   # view the data (should be comment to shorten run-time)

marker_name = re.findall(r'USDA_SNP\d{4}', data)  # a.Marker name
print(marker_name)
marker_ID = re.findall(r'[i]\d{5}\w{2}', data)  # b.Marker ID
marker_alleles_all = re.findall(r's=\w{1}\S\w{1}', data)  # both allels

# need to separate each allele from the above string
all_marker_len = len(marker_alleles_all)

print("number of different markers = ", all_marker_len)
print("typ: ", type(marker_alleles_all))

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


def find_snp_ind(text):  # return the index of the snp in the marker
    for ch in ['Y', 'M', 'K', 'R', 'W', 'S']:
        snp_ind = text.find(ch)
        if snp_ind != -1:
            break  # check if one of the 'ch' above is inside the string
    return snp_ind


for index in range(0, all_marker_len):  # finding the items
    Allele = marker_alleles_all[index]
    Allele_1 = Allele[2:3]   # c. Marker allele 1
    Allele_2 = Allele[4:5]   # c. Marker allele 2
    full_marker_seq_single_marker = full_marker_seq_raw[index]
    full_marker_seq = full_marker_seq_single_marker[2:]
    full_marker_seq_allel_1 = put_allel(full_marker_seq, Allele_1)  # e. Full marker sequence with allele 1
    full_marker_seq_allel_2 = put_allel(full_marker_seq, Allele_2)  # f. Full marker sequence with allele 2
    marker_length = len(full_marker_seq_allel_2)  # g. Marker length
    ind_snp = find_snp_ind(full_marker_seq)  # h. Position of SNP in the marker (0-based index)
    print(ind_snp)
    print(full_marker_seq_allel_1)
    print(full_marker_seq_allel_2)


# writing to file

dic = {"John": "john@example.com", "Mary": "mary@example.com"}  # dictionary

download_dir = "Markers_Data.csv"  # where you want the file to be downloaded to
csv = open(download_dir, "w")
columnTitleRow = "Marker_name, Marker_ID, Marker_Allele_1, Marker_Allele_2, Full_marker_sequence_Allele_1, \
 Full_marker_sequence_allele_2,Marker_length, Position_of_SNP\n"
csv.write(columnTitleRow)
'''
for key in dic.keys():
	name = key
	email = dic[key]
	row = name + "," + email + "\n"
	csv.write(row)


with open('Markers_Data.csv', 'wb') as csvfile:
    filewriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    filewriter.writerow([' Marker_name', 'Marker_ID', 'Marker_Allele_1', 'Marker_Allele_2', 'Full_marker_sequence_Allele_1'
,'Full_marker_sequence_allele_2','Marker_length', 'Position_of_SNP'])
    filewriter.writerow(['Derek', 'Software Developer'])
'''