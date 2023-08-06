import argparse
import json
from enum import Enum

class TaxaFilter(Enum):
    Kingdom = 0
    Phylum = 1
    Class = 2
    Genus = 3
    Species = 4
    Bacteria = 5
    Bacgroup = 6
    Order = 8


class HabitatFilter(Enum):
    Node = 0
    Bacteria = 1
    Kingdom = 2
    Pylum = 3
    Class = 4
    Genus = 5
    Species = 6
    Bacgroup = 7
    Taxid = 9
    Habitat = 10


habitat_filters_list = []
taxa_filters_list = []

class activity(Enum):
    B = "INFORMATION STORAGE AND PROCESSING Chromatin structure and dynamics"
    C = "METABOLISM Energy production and conversion"
    E = "METABOLISM Amino acid transport and metabolism"
    F = "METABOLISM Nucleotide transport and metabolism"
    G = "METABOLISM Carbohydrate transport and metabolism"
    H = "METABOLISM Coenzyme transport and metabolism"
    I = "METABOLISM Lipid transport and metabolism"
    J = "INFORMATION STORAGE AND PROCESSING Translation, ribosomal structure and biogenesis"
    K = "INFORMATION STORAGE AND PROCESSING Transcription"
    L = "INFORMATION STORAGE AND PROCESSING Replication, recombination and repair"
    O = "CELLULAR PROCESSES AND SIGNALING Posttranslational modification, protein turnover, chaperones"
    P = "METABOLISM Inorganic ion transport and metabolism"
    Q = "METABOLISM Secondary metabolites biosynthesis, transport and catabolism"
    R = "POORLY CHARACTERIZED General function prediction only"
    S = "POORLY CHARACTERIZED Function unknown"
    U = "CELLULAR PROCESSES AND SIGNALING Intracellular trafficking, secretion, and vesicular transport"
    V = "CELLULAR PROCESSES AND SIGNALING Defense mechanisms"
    W = "CELLULAR PROCESSES AND SIGNALING Extracellular structures"
    X = "MOBILOME Mobilome: prophages, transposons"
    Z = "CELLULAR PROCESSES AND SIGNALING Cytoskeleton"



def parse_args():
    """This function parses over the script parameters and extracts them."""
    parser = argparse.ArgumentParser(description='This parses over script arguments')
    parser.add_argument('-q', dest='genome_occurrences', required=True,
                        help='A number which defines the number of occurrences of cogx in genomes')
    parser.add_argument('-d', dest='motif_length', required=True, help='Length of motif')
    parser.add_argument('-x', dest='cogx', required=True, help='COGX to find')
    parser.add_argument('-pf', dest='plasmid_file', required=False, help='Plasmid file')
    parser.add_argument('-bf', dest='bacteria_file', required=False, help='Bacteria file')

    return parser.parse_args()

def check_taxa_filter(line_splitted_raw):
    found_strain_id = 0
    if taxa_filters_list:
        strain_id = int(line_splitted_raw[0].split('#')[-1])
        with open('taxa.txt', 'r') as taxa_file:
            taxa_file_lines = taxa_file.readlines()
            taxa_file_lines = taxa_file_lines[1:]
            for line in taxa_file_lines:
                if found_strain_id == 1:
                    break
                taxa_line_splitted_raw = line.strip('\n').split(',')
                if int(taxa_line_splitted_raw[6]) == strain_id:
                    found_strain_id = 1
                    for taxa_filter in taxa_filters_list:
                        if taxa_line_splitted_raw[taxa_filter[0]] != taxa_filter[1]:
                            return -1
    return 0


def check_habitat_filter(line_splitted_raw):
    found_bacteria_name = 0
    if habitat_filters_list:
        bacteria_name = line_splitted_raw[0].split('#')[3]
        with open('bactTaxa_Habitat.txt', 'r') as habitat_file:
            habitat_file_lines = habitat_file.readlines()
            habitat_file_lines = habitat_file_lines[1:]
            for line in habitat_file_lines:
                if found_bacteria_name == 1:
                    break
                habitat_line_splitted_raw = line.strip('\n').split(';')
                if habitat_line_splitted_raw[1] == bacteria_name:
                    found_bacteria_name = 1
                    for habitat_filter in habitat_filters_list:
                        if habitat_line_splitted_raw[habitat_filter[0]] != habitat_filter[1]:
                            return -1
    return 0


def find_cogx(cogx, d, q, plasmid_file, bacteria_file):
    """This function recieves:
    cogx- a specific cog to find inside the genomes.
    d - the minimnum length of a motif
    q - the minimum number of genomes in which the motif needs to be found.
    plasmid_file - the genome file containing words.
    bacteria_file - another genome file containing words.
    The function searches through the files line by line and searches for motifs
    with a minimum length of d that include cogx and the motif is found in more than q different genomes."""
    sorted_dmer_dict = {}
    dmer_dict = {}
    files = [plasmid_file, bacteria_file]
    for file in files:
        with open(file, 'r') as p_file:
            p_file_lines = p_file.readlines()
            # The final data structure which will be dumped to the json file.
            for line in p_file_lines:
                start_genome_index = line.index('_', 0) + 1
                genome_number = line[start_genome_index:start_genome_index + 6]  # Get genome number from line
                line_splitted_raw = line.split('\t')
                word = line_splitted_raw[1:-1]  # Extract only the cogs from the line
                word_length = len(word)  # Get word length
                if word_length < d:  # if smaller than d no need to look for motif, continue to next line
                    continue
                else:
                    for i in range(len(word) - d + 1):  # Iterate through word checking each motif available
                        start = i
                        end = i + d
                        try:
                            word.index(cogx, start, end)  # Check if cogx inside the motif, if not an exception is raised
                            if 'X' not in word[start:end]:  # make sure 'x not in motif, if so not relevant
                                motif = word[start:end]
                                motif_string = ' '.join(motif)
                                if check_taxa_filter(line_splitted_raw) == -1 or check_habitat_filter( #check filters
                                        line_splitted_raw) == -1:
                                    break
                                if motif_string in dmer_dict:  # motif was already found, add the genome number to the
                                    # dictionary, motif as key and genome number as value inside set
                                    dmer_dict[motif_string].add(genome_number)
                                else:  # First time motif is found, add to dictionary
                                    dmer_dict[motif_string] = set([genome_number])
                        except ValueError:  # cogx is not inside the 'kmer' inside the word, check next kmer inside word
                            continue
    for item in dmer_dict.items():  # Iterate over dictionary and add to final dictionary only the motifs that
        # appear in more than q genomes.
        motif_occurences_num = len(item[1])
        if motif_occurences_num >= q:
            if motif_occurences_num in sorted_dmer_dict:
                sorted_dmer_dict[motif_occurences_num].append(item[0])
            else:
                sorted_dmer_dict[motif_occurences_num] = [item[0]]
    return sorted_dmer_dict

def taxa_filters_func():
    numbers_chosen = []
    while input("Do you want to add a filter from the Taxa file to the COG search? y/n ").lower() == 'y':
        print("0: Kingdom\n"
              "1: Phylum\n"
              "2: Class\n"
              "3: Genus\n"
              "4: Species\n"
              "5: Bacteria\n"
              "6: Bacgroup\n"
              "8: Order\n\n")
        filter_num = int(input("Please write the number of the filter requested\n"))
        if filter_num in numbers_chosen:
            print("This filter has already been chosen\n")
            continue
        if filter_num not in [member.value for member in TaxaFilter]:
            print("Invalid number written")
            continue
        numbers_chosen.append(filter_num)
        filter_name = input(
            f"enter the value of the filter wanted for {TaxaFilter(filter_num).name}: ").lower().capitalize()
        taxa_filters_list.append((filter_num, filter_name))


def habitat_filters_func():
    numbers_chosen = []
    while input("Do you want to add a filter from the Habitat file to the COG search? y/n ").lower() == 'y':
        print("0: Node\n"
              "1: Bacteria\n"
              "2: Kingdom\n"
              "3: Pylum\n"
              "4: Class\n"
              "5: Genus\n"
              "6: Species\n"
              "7: Bacgroup\n"
              "9: Taxid\n"
              "10: Habitat\n\n")
        filter_num = int(input("Please write the number of the filter requested\n"))
        if filter_num in numbers_chosen:
            print("This filter has already been chosen\n")
            continue
        if filter_num not in [member.value for member in HabitatFilter]:
            print("Invalid number written")
            continue
        if filter_num not in [0, 8]:
            filter_name = input(
                f"enter the value of the filter wanted for {HabitatFilter(filter_num).name}: ").lower().capitalize()
        else:
            filter_name = int(input(
                f"enter the value number of the filter wanted for {HabitatFilter(filter_num).name}: "))
        habitat_filters_list.append((filter_num, filter_name))


def add_activity_to_dict(a_file, cog, a_dict, key):
    lines = a_file.readlines()
    for line in lines:
        if(cog == line[3:7]):
            if(line[8:9] not in a_dict.keys()):
                a_dict[line[8:9]] = 0   
                a_dict[line[8:9]] += key
                a_file.close
                break



def find_suspected_activity(d,cogx,dmers_dict):
    a_dict = {}
    cogs_num = 0
    for key in dmers_dict: #going through all the cogs found and characterizing their activity
        cogs_num += key * len(dmers_dict[key]) * (int(d) - 1) #calculating the number of cogs in the found words
        for motif in dmers_dict[key]:
            cogs = motif.split() 
            for cog in cogs:
                if cog == cogx:
                    continue
                else:
                    with open ("COG_INFO_TABLE.txt",'r') as a_file:
                        add_activity_to_dict(a_file, cog, a_dict, )
    sorted_dict = dict(sorted(a_dict.items(), key=lambda x: x[1],reverse=True))
    sorted_list= list(sorted_dict.items())
    final_dict={}
    print("Result Summary: after going through the adjacent cogs we suspect that those might be the cog's activities\n")

    for i in range(3):
        if(i<len(sorted_list)):
            precentage = round((sorted_list[i][1]/cogs_num)*100,1)
            print("The number",i,"most common activity found is:", activity[sorted_list[i][0]].value, "\nwith",precentage,"% hits\n")
            final_dict[activity[sorted_list[i][0]].value]= str(precentage)+"%"
    return final_dict
    




    
def main():
    q = input("please enter the minimal number of occurences of the cog in different genomes (q)")
    d = input("please enter the length of the word (d)")
    cogx = input ("please enter the cog you want to query (cogx)")
    taxa_filters_func()
    habitat_filters_func()
    print("\n"*3)
    dmers_dict=find_cogx(cogx, int(d), int(q), "cog_words_plasmid.txt","cog_words_bac.txt")
    sorted_dict = dict(  # sort the final dictionary by the q value (from small to big)
        sorted(dmers_dict.items(), reverse=True))
    activity_list=find_suspected_activity(d,cogx,dmers_dict)
    sorted_dict['Key defines dmer occurences in different genomes'] = 'Value is the dmer'
    sorted_cogs_json = json.dumps(sorted_dict, indent=4)
    activity_json = (json.dumps(activity_list, indent=4))
    with open(f'final_cogs_output_length{d}.json', 'w') as final_file:
        final_file.write(sorted_cogs_json)  # write final report4
    with open(f'final_suspected_{cogx}_activities_output.json', 'w') as final_file:
        final_file.write(activity_json)  # write final report4
    



main()

  