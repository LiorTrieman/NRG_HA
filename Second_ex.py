# ---------Part 2 - genotyping data------------- #

""" The TAMU array was used to genotype hundreds of cotton varieties. Genotyping results are
provided in the file GT_AH_DiversityAnalysisDataforCottonGen.xlsx, which can be
downloaded from:
https://nrgene-candidatetask.s3.amazonaws.com/GT_AH_DiversityAnalysisDataforCottonGen.xlsx
In this section, we’ll parse this file and extract some statistics. You may convert the xls into
a csv file for easier parsing.
1. Write a script that will read the genotyping file and outputs a file with commaseparated table. Each row of the table should represent one sample and have the
following information:
a. Sample name - e.g. AH-097
b. Sample ID (entry) - e.g. 392
c. Sample description (designation) - e.g. G. amourianum
d. Fraction of markers successfully genotyped - e.g. 0.92
e. Fraction of heterozygous genotypes (out of the successfully genotyped markers) -
e.g. 0.11
Ensure that genotypes match the possible alleles per marker, as defined in Part 1. If
you find mismatches, discard the genotype and issue an appropriate warning."""
import pandas as pd
import re  # in order to use regular expressions
import numpy as np  # arrange data
import csv
import matplotlib.pyplot as plt  # create plots and graphs
from collections import Counter


df = pd.read_csv('GT_AH_CSV', delimiter=',')  # read csv as df

# get the lists from the df:

sample_name_list = list(df.columns.values)  # headers
df_sample_ID_list = df.iloc[0, :].copy()
df_sample_description_list = df.iloc[1, :].copy()
df_codons_list = df.iloc[:, 1].copy()
print(df_codons_list.head(11))
# print(sample_name_list)
sample_ID_list = df_sample_ID_list.values.tolist()
# print(sample_ID_list)
sample_description_list = df_sample_description_list.values.tolist()
# print(sample_description_list)
codons_list = df_codons_list.values.tolist()
codons_list_no_head = codons_list[2:]
#print(codons_list_no_head)
#print(len(codons_list_no_head))
snp_objects = Counter(codons_list_no_head).keys()  # equals to list(set(words))
K_snp_counts = codons_list_no_head.count('K')  # counts the elements' appearances
M_snp_counts = codons_list_no_head.count('M')  # counts the elements' appearances
W_snp_counts = codons_list_no_head.count('W')  # counts the elements' appearances
R_snp_counts = codons_list_no_head.count('R')  # counts the elements' appearances
Y_snp_counts = codons_list_no_head.count('Y')  # counts the elements' appearances
T_codon_counts = codons_list_no_head.count('T')  # counts the elements' appearances
G_codon_counts = codons_list_no_head.count('G')  # counts the elements' appearances
C_codon_counts = codons_list_no_head.count('C')  # counts the elements' appearances
A_codon_counts = codons_list_no_head.count('A')  # counts the elements' appearances
total_codons = A_codon_counts + C_codon_counts + G_codon_counts+ T_codon_counts
total_snps = Y_snp_counts + R_snp_counts + W_snp_counts + M_snp_counts + K_snp_counts

Fraction_success_markers = total_codons / (total_codons + total_snps)
print(Fraction_success_markers)


def get_counter_snp(codons_list_no_head):
    K_snp_counts = codons_list_no_head.count('K')  # counts the elements' appearances
    M_snp_counts = codons_list_no_head.count('M')  # counts the elements' appearances
    W_snp_counts = codons_list_no_head.count('W')  # counts the elements' appearances
    R_snp_counts = codons_list_no_head.count('R')  # counts the elements' appearances
    Y_snp_counts = codons_list_no_head.count('Y')  # counts the elements' appearances
    total_snps = Y_snp_counts + R_snp_counts + W_snp_counts + M_snp_counts + K_snp_counts
    return total_snps


def get_counter_codons(codons_list_no_head):
    T_codon_counts = codons_list_no_head.count('T')  # counts the elements' appearances
    G_codon_counts = codons_list_no_head.count('G')  # counts the elements' appearances
    C_codon_counts = codons_list_no_head.count('C')  # counts the elements' appearances
    A_codon_counts = codons_list_no_head.count('A')  # counts the elements' appearances
    total_codons = A_codon_counts + C_codon_counts + G_codon_counts + T_codon_counts
    return total_codons


def write_raw_to_csv(csv, index, Fraction_success_markers):
    row = sample_name_list[index] + "," + sample_ID_list[index] + "," + sample_description_list[index] + "," \
          + str(Fraction_success_markers) + "," + "\n"
    csv.write(row)  # not the fastest way of doing it..



download_dir = "Part_2_Samples_Data.csv"  # where you want the file to be downloaded to
csv = open(download_dir, "w")
columnTitleRow = "Sample_Name, Sample_ID, Sample_Description, Fraction_successful, Fraction_genotype, \n"
csv.write(columnTitleRow)

# now apply on all columns of the df
Fract_success_markers_list = [] # init fraction list
for index in range(1, len(sample_name_list)): #starting with 10 samples
    df_codons_list = df.iloc[:, index].copy()
    codons_list = df_codons_list.values.tolist()
    codons_list_no_head = codons_list[2:]
    total_codons = get_counter_codons(codons_list_no_head)
    total_snps = get_counter_snp(codons_list_no_head)
    Fraction_success_markers = total_codons / (total_codons + total_snps)
    Fract_success_markers_list += str(Fraction_success_markers)  # convert values to str and get the list of fractions
    write_raw_to_csv(csv, index, Fraction_success_markers)
# WRITE THE EXTRACTED DATA TO A CSV FILE
# Creating a CSV file


"""
for index in range(0, len(sample_name_list)):  # finding the items
    row = sample_name_list[index] + "," + sample_ID_list[index] + "," + sample_description_list[index] + "," \
          + str(Fraction_success_markers) + "," + "\n"
    csv.write(row)  # not the fastest way of doing it..
"""

