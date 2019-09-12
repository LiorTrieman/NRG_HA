# PART 2 - GENOTYPING DATA  #
# ------------------------- #
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
from collections import Counter
import matplotlib.pyplot as plt  # create plots and graphs


def get_marker_list(data_frame):  # get the lists from the df:
    df_marker_list = data_frame.iloc[:, 0].copy()
    marker_list = df_marker_list.values.tolist()
    marker_list_no_head = marker_list[2:]
    return marker_list_no_head


def get_counter_codons(seq_list_no_head):
    T_codon_counts = seq_list_no_head.count('T')  # counts the elements' appearances
    G_codon_counts = seq_list_no_head.count('G')  # counts the elements' appearances
    C_codon_counts = seq_list_no_head.count('C')  # counts the elements' appearances
    A_codon_counts = seq_list_no_head.count('A')  # counts the elements' appearances
    total_codons = A_codon_counts + C_codon_counts + G_codon_counts + T_codon_counts
    return total_codons


def get_counter_failed(seq_list_no_head):
    failed_counts = seq_list_no_head.count('failed')  # counts the failed appearances
    return failed_counts


def write_raw_to_csv(csv, index, Fraction_success_markers, fraction_hetro_markers):
    row = sample_name_list[index] + "," + sample_ID_list[index] + "," + sample_description_list[index] + "," \
          + str(Fraction_success_markers) + "," + str(fraction_hetro_markers) + "\n"
    csv.write(row)


def get_counter_snp(seq_list_no_head):
    K_snp_counts = seq_list_no_head.count('K')  # counts the elements' appearances
    M_snp_counts = seq_list_no_head.count('M')  # counts the elements' appearances
    W_snp_counts = seq_list_no_head.count('W')  # counts the elements' appearances
    R_snp_counts = seq_list_no_head.count('R')  # counts the elements' appearances
    Y_snp_counts = seq_list_no_head.count('Y')  # counts the elements' appearances
    total_snps = Y_snp_counts + R_snp_counts + W_snp_counts + M_snp_counts + K_snp_counts
    return total_snps

#   QUESTIONS 1 - CREATE A CSV FILE FROM GENOTYPING DATA  #
# ------------------------------------------------------- #


if __name__ == "__main__":

    with open('TAMU_SNP63K_69997.fasta', 'r') as file:  # read the fasta file
        data = file.read().replace('\n', '')

    marker_alleles_all = re.findall(r's=\w{1}\S\w{1}', data)  # both alleles
    # need to separate each allele from the above string
    all_marker_len = len(marker_alleles_all)
    Allele_1_list = []
    Allele_2_list = []
    for index in range(0, all_marker_len):  # finding the items
        Allele = marker_alleles_all[index]
        Allele_1 = Allele[2:3]  # c. Marker allele 1
        Allele_2 = Allele[4:5]  # d. Marker allele 2
        Allele_1_list += Allele_1  # use to check fraction of genotype
        Allele_2_list += Allele_2  # use to check fraction of genotype

    df = pd.read_csv('GT_AH_CSV', delimiter=',')  # read csv as df

    get_marker_list(df) # get list of markers to check for hetrozigoty
    sample_name_list = list(df.columns.values)  # headers
    df_sample_ID_list = df.iloc[0, :].copy()
    df_sample_description_list = df.iloc[1, :].copy()
    df_seq_list = df.iloc[:, 1].copy()
    sample_ID_list = df_sample_ID_list.values.tolist()
    sample_ID_list_keys = Counter(sample_ID_list).keys()  # equals to list(set(words))
    sample_ID_list_values = Counter(sample_ID_list).values()  # counts the elements' appearances

    sample_description_list = df_sample_description_list.values.tolist()
    seq_list = df_seq_list.values.tolist()
    seq_list_no_head = seq_list[2:]

    snp_objects = Counter(seq_list_no_head).keys()  # equals to list(set(words))
    K_snp_counts = seq_list_no_head.count('K')  # counts the elements' appearances
    M_snp_counts = seq_list_no_head.count('M')  # counts the elements' appearances
    W_snp_counts = seq_list_no_head.count('W')  # counts the elements' appearances
    R_snp_counts = seq_list_no_head.count('R')  # counts the elements' appearances
    Y_snp_counts = seq_list_no_head.count('Y')  # counts the elements' appearances
    T_codon_counts = seq_list_no_head.count('T')  # counts the elements' appearances
    G_codon_counts = seq_list_no_head.count('G')  # counts the elements' appearances
    C_codon_counts = seq_list_no_head.count('C')  # counts the elements' appearances
    A_codon_counts = seq_list_no_head.count('A')  # counts the elements' appearances
    total_codons = A_codon_counts + C_codon_counts + G_codon_counts + T_codon_counts
    total_snps = Y_snp_counts + R_snp_counts + W_snp_counts + M_snp_counts + K_snp_counts

    Fraction_success_markers = total_codons / (total_codons + total_snps)
    Fraction_hetro_markers = total_snps / (total_codons + total_snps)

    download_dir = "Part_2_Sample_Data.csv"  # where you want the file to be downloaded to
    csv = open(download_dir, "w")
    columnTitleRow = "Sample_Name, Sample_ID, Sample_Description, Fraction_successful, Fraction_genotype, \n"
    csv.write(columnTitleRow)

    # now apply on all columns of the df
    fract_success_markers_list = []  # init fraction list
    fraction_hetro_markers_list = []  # init fraction list
    for index in range(1, len(sample_name_list)):  # starting with 10 samples
        df_seq_list = df.iloc[:, index].copy()
        seq_list = df_seq_list.values.tolist()
        seq_list_no_head = seq_list[2:]
        total_codons = get_counter_codons(seq_list_no_head)
        total_snps = get_counter_snp(seq_list_no_head)
        fraction_success_markers = total_codons / (total_codons + total_snps)
        fraction_hetro_markers = total_snps / (total_codons + total_snps)
        fract_success_markers_list += str(fraction_success_markers)  # convert values to str,get the list of fractions
        fraction_hetro_markers_list += str(fraction_hetro_markers)  # convert values to str,get the list of fractions
        write_raw_to_csv(csv, index, fraction_success_markers, fraction_hetro_markers)

    #   QUESTIONS 2 - PLOT HISTOGRAMS OF THE FRACTION OF MISSING DATA BY SAMPLE AND BY MARKER  #
    # ---------------------------------------------------------------------------------------- #

    # GET MISSING DATA FRACTION OUT OF EACH SAMPLE
    total_failed_sample_list = []
    fraction_failed_per_sample_list = []
    number_of_markers = all_marker_len
    print("number_of_markers: ",number_of_markers)
    for index in range(1, len(sample_name_list)):  # starting with 10 samples #
        df_seq_list = df.iloc[:, index].copy()
        seq_list = df_seq_list.values.tolist()
        seq_list_no_head = seq_list[2:]
        total_failed_per_sample = get_counter_failed(seq_list_no_head)
        total_failed_sample_list.append(total_failed_per_sample)
        fraction_failed_per_sample = total_failed_per_sample/number_of_markers
        fraction_failed_per_sample_list.append(fraction_failed_per_sample)

    # GET MISSING DATA FRACTION OUT OF EACH MARKER
    total_failed_marker_list = []
    fraction_failed_per_marker_list = []
    number_of_samples = len(sample_name_list)
    marker_list_no_head = get_marker_list(df)
    for index in range(1,len(marker_list_no_head)):  #
        df_seq_list = df.iloc[index, :].copy()
        seq_list = df_seq_list.values.tolist()
        seq_list_no_head = seq_list[1:]
        total_failed_per_marker = get_counter_failed(seq_list_no_head)
        total_failed_marker_list.append(total_failed_per_marker)
        fraction_failed_per_marker = total_failed_per_marker/number_of_samples
        fraction_failed_per_marker_list.append(fraction_failed_per_marker)
# PLOTTING THE COUNTS
    failed_keys = Counter(total_failed_marker_list).keys()  # equals to list(set(words))
    failed_counts = Counter(total_failed_marker_list).values()  # counts the elements' appearances
    plt.plot(failed_keys, failed_counts, "*")
    plt.title("Fraction of failed per marker")
    plt.show()
    failed_keys = Counter(fraction_failed_per_sample_list).keys()  # equals to list(set(words))
    failed_counts = Counter(fraction_failed_per_sample_list).values()  # counts the elements' appearances
    plt.plot(failed_keys, failed_counts, "*")
    plt.title("Fraction of failed per sample")
    plt.show()

    plt.hist(fraction_failed_per_marker_list, density=False, histtype='bar', align='mid', orientation='vertical')
    plt.title('Fraction of Failed  Per Marker')
    plt.show()
    plt.figure(2)
    plt.hist(fraction_failed_per_sample_list, density=False, histtype='bar', align='mid', orientation='vertical')
    plt.title('Fraction of Failed  Per Sample')
    plt.show()






