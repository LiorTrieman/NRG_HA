# QUESTION #1#
#
#
#
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
import numpy as np # arrange data
import matplotlib.pyplot as plt  # create plots and graphs
from collections import Counter

with open('TAMU_SNP63K_69997.fasta', 'r') as file:  # read the fasta file
    data = file.read().replace('\n', '')
# print(type(data))  # view data type

marker_name = re.findall(r'name=\S{1,50}', data)  # a.Marker name
# print("len marker name: ", len(marker_name))  # to verify that is finds all names
marker_ID = re.findall(r'63K_i\d{5}\w{2}', data)  # b.Marker ID
# print("len marker id: ", len(marker_ID))

marker_alleles_all = re.findall(r's=\w{1}\S\w{1}', data)  # both alleles
# print(len(marker_alleles_all))
# need to separate each allele from the above string
all_marker_len = len(marker_alleles_all)

print("number of different markers = ", all_marker_len)
print("typ: ", type(marker_alleles_all))

# to find  Full marker sequence we should search for :spaceX2, {A,G,C,T, and Y,M,K,R,W,S} and stop at "enter" or "<"

full_marker_seq_raw = re.findall(r'\s\s[A, G, C, T, Y, M, K, R, W, S]*', data)  # Full marker sequence with allele 1
# print(full_marker_seq_raw[0:10])
full_marker_seq_raw_len = len(full_marker_seq_raw)
# print("len: ", full_marker_seq_raw_len)

# Creating a CSV file
download_dir = "Markers_Data.csv"  # where you want the file to be downloaded to
csv = open(download_dir, "w")
columnTitleRow = "Marker_name, Marker_ID, Marker_Allele_1, Marker_Allele_2, Full_marker_sequence_Allele_1, \
 Full_marker_sequence_allele_2,Marker_length, Position_of_SNP\n"
csv.write(columnTitleRow)

flag_more_than_one_snp = [0] * len(marker_alleles_all)


def put_allele(text, allele):  # replace one of [Y, M, K, R, W, S] with the alleles
    for ch in ['Y', 'M', 'K', 'R', 'W', 'S']:
        if ch in text:
            text = text.replace(ch, allele)  # replace with one of the alleles
    return text


def find_snp_ind(text):  # return the index of the snp in the marker
    for ch in ['Y', 'M', 'K', 'R', 'W', 'S']:
        snp_ind = text.find(ch)
        if snp_ind != -1:
            break  # check if one of the 'ch' above is inside the string
    return snp_ind
# HAVEN'T FOUND A MARKER WITH MORE THAN ONE WILDCARD CHARACTER


# print(all_marker_len)
for index in range(0, all_marker_len):  # finding the items
    marker_name_current_raw = marker_name[index]
    marker_name_current = marker_name_current_raw[5:]
    marker_ID_current_raw = marker_ID[index]
    marker_ID_current = marker_ID_current_raw[4:]
    Allele = marker_alleles_all[index]
    Allele_1 = Allele[2:3]   # c. Marker allele 1
    Allele_2 = Allele[4:5]   # c. Marker allele 2
    full_marker_seq_single_marker = full_marker_seq_raw[index]
    full_marker_seq = full_marker_seq_single_marker[2:]
    full_marker_seq_allele_1 = put_allele(full_marker_seq, Allele_1)  # e. Full marker sequence with allele 1
    full_marker_seq_allele_2 = put_allele(full_marker_seq, Allele_2)  # f. Full marker sequence with allele 2
    marker_length = len(full_marker_seq_allele_2)  # g. Marker length
    ind_snp = find_snp_ind(full_marker_seq)  # h. Position of SNP in the marker (0-based index)
    # flag the seq with more than one snp:
    ''' if len(str(ind_snp)) > 1:
        flag_more_than_one_snp[index] = 1
    else:
        flag_more_than_one_snp[index] = 0 '''
    row = marker_name_current + "," + marker_ID_current + "," + Allele_1 + "," + Allele_2 + "," + full_marker_seq_allele_1 \
          + "," + full_marker_seq_allele_2 + "," + str(marker_length) + "," + str( ind_snp) + "," + "\n"
    csv.write(row) # not the fastest way of doing it..

# print(flag_more_than_one_snp)

# PART #1 QUESTION #2
#
#
''' Use the output of the above script to simplify the fasta - write a script that will read
the output and print two fasta files with full marker sequences, one file for allele 1 and the
other for allele 2. The fasta headers should only include the marker name#
'''

# first create a list of marker_names and a list of marker_seq for each allele
# create a dict from these two lists
# write then into fasta file
list_name = []
list_seq = []

for index in range(0, all_marker_len):  # finding the items
    marker_name_current_raw = marker_name[index]
    marker_name_current = marker_name_current_raw[5:]
    list_name.append(marker_name_current)
print(list_name)
print(len(list_name))


o_file = open("fasta_allele_1.txt", "w")
# write a title of names:
o_file.write(">")
for index in range(len(list_name)):
    o_file.write(list_name[index] + ",")
# write allele 1 seq
for index in range(len(list_name)):
    Allele = marker_alleles_all[index]
    Allele_1 = Allele[2:3]  # c. Marker allele 1
    full_marker_seq_single_marker = full_marker_seq_raw[index]
    full_marker_seq = full_marker_seq_single_marker[2:]
    full_marker_seq_allele_1 = put_allele(full_marker_seq, Allele_1)  # e. Full marker sequence with allele 1
    o_file.write(">" + full_marker_seq_allele_1 + "\n")
o_file.close()
o_file = open("fasta_allele_2.txt", "w")
# write a title of names:
o_file.write(">")
for index in range(len(list_name)):
    o_file.write(list_name[index] + ",")
# write allele 2 seq
for index in range(len(list_name)):
    Allele = marker_alleles_all[index]
    Allele_2 = Allele[4:5]  # c. Marker allele 2
    full_marker_seq_single_marker = full_marker_seq_raw[index]
    full_marker_seq = full_marker_seq_single_marker[2:]
    full_marker_seq_allele_2 = put_allele(full_marker_seq, Allele_2)  # f. Full marker sequence with allele 2
    o_file.write(">" + full_marker_seq_allele_2 + "\n")
o_file.close()


#  PLOTING THE STATISTICS OF THE DATA

# Bar Chart for Frequencies of SNP types (C/T, C/G. C/A, A/T, A/G, G/T)
# marker_alleles_all -raw list of all SNP types, need to extract freqs
num_total_snp = len(marker_alleles_all)
list_of_snp = []


def cut_head_of_snp_string(raw_SNP):
    SNP = raw_SNP[2:]
    return SNP


for index in range(0, num_total_snp):
    snp_type = cut_head_of_snp_string(marker_alleles_all[index])
    list_of_snp.append(snp_type)

snp_freqs = []
snp_objects = Counter(list_of_snp).keys()  # equals to list(set(words))
snp_counts = Counter(list_of_snp).values()  # counts the elements' appearances
for index in snp_counts:
    freq = index/num_total_snp
    snp_freqs.append(freq)
y_pos = np.arange(len(snp_objects))
plt.figure(1)  # new figure
fig = plt.bar(y_pos, snp_freqs, align='center', alpha=0.5)
plt.xticks(y_pos, snp_objects)
plt.ylabel('Frequencies')
plt.title('Frequencies of SNP types')
plt.show()
# Marker length - Histogram
list_of_marker_length = []


def cut_head_of_marker_seq_string(raw_str):
    length = raw_str[2:]
    return length


for index in range(0, num_total_snp):
    marker_str = cut_head_of_marker_seq_string(full_marker_seq_raw[index])
    marker_len = len(marker_str)
    list_of_marker_length.append(marker_len)

# print(list_of_marker_length)  # test print
# print(len(list_of_marker_length))  # test print
plt.figure(2)
fig_lengths = plt.hist(list_of_marker_length, range=(0, 400), density=False, histtype='bar', align='mid', orientation='vertical')
plt.ylabel('Marker Lengths counts')
plt.xlabel('Marker Lengths (number of Nucleotide)')
plt.title('Marker Lengths Histogram')
plt.show()
''' bins=None, , density=None, weights=None, cumulative=False,\
                       bottom=None, , rwidth=None, log=False, color=None, label=None,\
                       stacked=False, normed=None, *, data=None, **kwargs)'''

# Position of SNP in marker sequence - normalized to 1 (pos/len)
# correlation between snp typp (c/t or c/g for example) and marker length.
# correlation between snp location and marker length


# ----------Part 2 - genotyping data-----------##
# ---------------------------------------------##


