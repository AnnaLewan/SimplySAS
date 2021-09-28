#!/usr/bin/env python3

import sys
import os
import pandas as pd
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import pairwise2
import zipfile
import tarfile
import gzip
from collections import defaultdict
import uuid
import subprocess
import shutil
import argparse
from multiprocessing import Pool
import time

parser = argparse.ArgumentParser(prog='SimplySAS', description="Reimplementation of AmpliSAS. Performs genotyping of amplicon sequencing data by clustering errors and filtering artefacts following a 3-step based analysis: de-multiplexing, clustering and filtering.")

parser.add_argument("-i",'--INPUT',type=str,required=True,help="Input set of FASTQ or FASTA files packed into a unique .ZIP or .TAR.GZ file.")
parser.add_argument("-m",'--METHOD',type=str,required=True,help="Method used to obtain reads (old or new). Old - demultiplexing needed, new - without demultiplexing.")
parser.add_argument("-bp",'--BARPRI',type=str,required=False,help="CSV file with primer and barcode data for old method.")
parser.add_argument("-p",'--PRIMER',type=str,required=False,help="CSV file with primer data for new method.")
parser.add_argument("-o",'--OUTPUT',type=str,required=True,help="Output folder name.")
parser.add_argument("-el", '--EXPLEN',type=int,required=False,help="Expected length of the marker sequence.")
parser.add_argument("-se", '--SUBERROR',type=float,required=False,help="Threshold for permissible substitution error (%).")
parser.add_argument("-df", '--DOMFREQ',type=float,required=False,help="Minimum sequence frequency respect to the dominant (%).")
parser.add_argument("-af", '--AMFREQ',type=float,required=False,help="Minimum sequence frequency per-amplicon (%).")
parser.add_argument("-ch", '--CHIMERA',type=int,required=False,help="Minimal length of match within sequences to consider as chimera (bp).")
parser.add_argument("-al", '--MAXALL',type=int,required=False,help="Maximal number of allels in one amplicon.")
parser.add_argument("-ad", '--MINDEPTH',type=int,required=False,help="Minimal depth of amplicon to not be discarded.")
parser.add_argument("-ma",'--MINALL',type=int,required=False,help="Minimal allel depth to not be discarded in filtering.")
parser.add_argument("-t",'--THREADS',type=int,required=False,help="Number of threads used in analysis.")

args = parser.parse_args()

### Checking parametres ###

# Required
if args.INPUT == None:
	raise Exception('No input file')
else:
	input_multifile = args.INPUT
		
if os.path.isfile(input_multifile) != True:
    raise Exception('Input multifile does not exist')
    
    
if args.METHOD == None:
	raise Exception('Method not selected')
else:
	method = args.METHOD


if args.BARPRI == None and args.PRIMER == None:
	raise Exception('Primer/barcode file does not exist')

    
if args.BARPRI == None:
	barcode_file = None
else:
	barcode_file = args.BARPRI
	
if args.PRIMER == None:
	primers_file = None
else:
	primers_file = args.PRIMER
    

if args.OUTPUT == None:
	raise Exception('No output directory')
else:
	output_dir = args.OUTPUT
    

# Clustering parametres

if args.SUBERROR == None:
    substitution_error_threshold = 1
else:
    substitution_error_threshold = args.SUBERROR #threshold dla dopuszczalnego błędu substytucji aby odczyt mógł zostać przypisany do danego klastra (w %)


if args.DOMFREQ == None:
    min_dominant_frequency_threshold = 25
else:
	min_dominant_frequency_threshold = args.DOMFREQ #granica częstości występowania w odniesieniu do sekwencji dominującej


if args.AMFREQ == None:
    min_amplicon_seq_frequency = 0.25
else:
	min_amplicon_seq_frequency = args.AMFREQ # granica częstości występowania w odniesieniu do wszystkich sekwencji w amplikonie

# Filtering parameters
if args.CHIMERA == None:
    min_chimera_length = 10
else:
	min_chimera_length = args.CHIMERA


if args.MAXALL == None:
    max_allele_number = 75
else:
	max_allele_number = args.MAXALL


if args.MINDEPTH == None:
    min_amplicon_depth = 100
else:
	min_amplicon_depth = args.MINDEPTH


if args.MINALL == None:
    min_allel_depth = 2
else:
    min_allel_depth = args.MINALL
    

# Threads parameter
if args.THREADS == None:
	threads = 4
else:
	threads = args.THREADS
	

### Folder creating ###

# Główny folder
p1 = f'{output_dir}/ampliSAS_analysis'
os.makedirs(p1)

# Folder wstępnej analizy
prework_dir = f'{p1}/prework'
os.mkdir(prework_dir)

# Folder do wypakowania plików
pre_input = f'{prework_dir}/pre_input'
os.mkdir(pre_input)

# Folder plików po merge
merged = f'{prework_dir}/merged'
os.mkdir(merged)

# Folder plików po przyczinaniu
trimmed = f'{prework_dir}/trimmed'
os.mkdir(trimmed)

# Folder plików po usunięciu starterów
clean = f'{prework_dir}/clean'
os.mkdir(clean)

# Folder plików po usunięciu metek i demultipleksowaniu
disbarred_dir = f'{clean}/disbarred'
os.mkdir(disbarred_dir)

# 
p7 = f'{p1}/fasta'
os.mkdir(p7)

p8 = f'{p1}/amplicons_filtered_seqs' # Tu będa znajdowały się foldery dla każdego amplikonu zawierające każdą sekwencję w amplikonie jako osobny plik fasta. Sekwencje te są już przefiltrowane względem długości, głebi itd
os.mkdir(p8)

p9 = f'{p1}/clusters'
os.mkdir(p9)

# Folder plików
p10 = f'{p1}/filtered'
os.mkdir(p10)


### Files unpacking ###

shutil.unpack_archive(input_multifile, pre_input)


###############################################
################# CLASSES #####################

class Prework:
    def __init__(self):
        pass
    
    
    def pear(self, file_left, file_right, name, outdir, threads):
        output = f'{outdir}/{name}'
        command = f'pear -f {file_left} -r {file_right} -o {output} -v 5 -j {threads}'
        process = subprocess.check_output(command, shell=True)
        
    
    def grep(self, file, name, outdir):        
        command = f'cat {file} | seqkit grep -i -s -v -p N > {outdir}/{name}_Nout.fastq'
        process = subprocess.check_output(command, shell=True)
        
    def cut_primers(self, file, name, forward, reverse, region_f, region_r, outdir):
        command = f'seqkit amplicon {file} -F {forward} -R {reverse} -r {region_f}:{region_r} -o {outdir}/{name}_trimmed.fastq'
        process = subprocess.check_output(command, shell=True)
        
        
        
class Amplicon:
    def __init__(self, number):
        self.number = number
        self.clusters = None
        self.amplicon_table = pd.DataFrame(columns = ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'])
        
            
    def parse_sequence_file(self, inputfastq, method_end):         

        dedup_records = defaultdict(list) #tworzy słownik, dla których wartościami jest lista

        for record in SeqIO.parse(inputfastq, method_end):
            dedup_records[str(record.seq)].append(record.id) #tworzy słownik, dla którego key=seq, a value=id. Na tym etapie dochodzi do odnalezienia duplikatów - jest ich tyle wartości id
        for seq, ids in sorted(dedup_records.items(), key=lambda t: len(t[1]), reverse=True): #ta linijka służy przesortowaniu sekwencji pod względem ich głębokości                
            id_hold = ids[0]
            id_len = len(ids)
            seq_len = len(seq)
            unique_hash = str(uuid.uuid1())
            self.amplicon_table = self.amplicon_table.append({'Hash' : unique_hash, 'ID' : id_hold, 'Sequence' : seq, 'Depth' : id_len, 'Length' : seq_len}, ignore_index = True)
        
        depth_of_amplicon = self.amplicon_table['Depth'].sum() #Wylicza głębie całego amplikonu
        
        self.amplicon_table['Frequency'] = self.amplicon_table['Depth']/depth_of_amplicon*100
                    
        return self.amplicon_table

    
    def find_minmax(self):
        if args.EXPLEN == None:
            expected_length = self.amplicon_table.loc[0, 'Length']
        else:
            expected_length = args.EXPLEN
        
        max_length = expected_length + 15
        min_length = expected_length - 15
        
        self.amplicon_table.drop(self.amplicon_table[self.amplicon_table.Length < min_length].index, inplace=True)  #odrzuca sekwencje o długości niższej niż zadana przez użytkownika
        self.amplicon_table.drop(self.amplicon_table[self.amplicon_table.Length > max_length].index, inplace=True)  #odrzuca sekwencje o długości wyższej niż zadana przez użytkownika
        self.amplicon_table.reset_index(drop=True, inplace=True)
        
    
    def create_fasta(self, inputs_dir):
        outpath = f'{inputs_dir}/{self.number}_amplicon.fasta'
        with open(outpath, 'w') as outfile:
            for index, row in self.amplicon_table[ ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
                outfile.write('>' + row.Hash + ' | ' + row.ID + ' | depth: ' + str(row.Depth) + ' | length: ' + str(row.Length) + ' | frequency per amplicon: ' + str(row.Frequency) + '\n' + row.Sequence + '\n')
                
                
    def create_seq_fasta(self, seqs_fasta_dir):
        for index, row in self.amplicon_table[ ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
            outpath = f'{seqs_fasta_dir}/{row.Hash}.fasta'
            with open(outpath, 'w') as outfile:
                outfile.write('>' + row.Hash + ',' + row.ID + ',' + str(row.Depth) + ',' + str(row.Length) + ',' + str(row.Frequency) + '\n' + row.Sequence + '\n')


class Cluster:
    def __init__(self, amplicon):
        self.dominant_seq = None
        self.amplicon = amplicon
                               
                
    def seq2seq(self, seq1, seq2, substitution_error_threshold, number):
        command = f'./fogsaa {seq1} {seq2} 1 0 1 -1 -1'
        process = subprocess.check_output(command, shell=True) #Przykładowy wynik: b'176\n176\nElapsed time: 1 milliseconds\ntotal nodes expanded==176\n\nscore= 174\n'
        process = process.decode('utf-8') #Zamiana bajtów na str

                
        score = int(re.sub('score=\ ', '', ' '.join(map(str, re.findall('score=\ .[0-9]*', process))))) # Wyodrębnienie score dla alignmentu
        #print(f'Score: {score}')
        len_seqs = int(re.sub('score=\ ', '', ' '.join(map(str, re.findall('(?<=[0-9]\s)[0-9]+', process))))) # Wyodrębnia długość porównywanych sekwencji - założenie, że porównywane sekwencje są tylko tej samej długości
        
        
        global substitution_error
        substitution_error = (len_seqs - score)/len_seqs*100 # Wylicza poziom błedów dla danego alignmentu
        
        
    def multiple_aline_seqs(self, cluster_df, clusters_dir, number):
        #Tworzenie wspólnej fasty dla klastra
        outpath = f'{clusters_dir}/{number}_cluster.fasta'
        with open(outpath, 'w') as outfile:
            for index, row in cluster_df[ ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
                outfile.write('>' + row.Hash + ' | ' + row.ID + ' | depth: ' + str(row.Depth) + ' | length: ' + str(row.Length) + ' | frequency per amplicon: ' + str(row.Frequency) + '\n' + row.Sequence + '\n')
                                
        # Stworzenie dataframe na sekwencję konsensusową z danymi
        global consensus_df
        
        consensus_df = pd.DataFrame(columns = ['ID', 'Sequence', 'Depth', 'Length'])
        id_cluster = f'cluster_{number}'
        len_cluster = cluster_df.iloc[0, 4]
        depth_cluster = cluster_df['Depth'].sum()
               
        #Ustalenie konsensusowej
        align = AlignIO.read(outpath, 'fasta')
        summary_align = AlignInfo.SummaryInfo(align)
        consensus = str(summary_align.dumb_consensus())
    
        consensus_df = consensus_df.append({'ID' : id_cluster, 'Depth' : depth_cluster, 'Length' : len_cluster, 'Sequence' : consensus}, ignore_index = True)



class Filter:
    def __init__(self):
        pass

        
    def is_chimera(self, cons_seqs_df):
        seqs_count = len(cons_seqs_df)
        cons_seqs_df['Chimera'] = 'Not'
        if seqs_count >= 3:
            reversed_df = cons_seqs_df.sort_values(['Depth'])
            
            indexer = [i for i in range(0, seqs_count)]
            #Do ustalenia
            identity_threshold = 30
            
            
            for j in indexer:
                for i in indexer:
                    seq1 = str(reversed_df.at[j, 'Sequence'])
                    seq2 = str(cons_seqs_df.at[i, 'Sequence'])
                    alignment = str(pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True))
                    search1 = re.search('(?<=seqB=\')[ACTGX\-]+\'', alignment)
                    search2 = re.search('(?<=score=)[0-9]+', alignment)
                    search3 = re.search('(?<=end=)[0-9]+', alignment)
                    aligned_seq = search1.group(0)
                    score = int(search2.group(0))
                    align_length = int(search3.group(0))
                    
                    #Wylicza identity
                    identity = score/align_length*100
                    if identity < identity_threshold:
                        seq2_list = list(aligned_seq)
                        # Tworzy binary score
                        binary_score = []
                        for base in seq2_list:
                            if base == 'A' or base == 'T' or base == 'C' or base == 'G':
                                binary_score.append(1)
                            else:
                                binary_score.append(0)
                                
                        binary_score_right_end = binary_score[:min_chimera_length] #Wycinek długości chimera_length bp z prawej
                        binary_score_left_end = binary_score[-min_chimera_length:] #Wycinek długości chimera_length bp z lewej
                        
                        #Sprawdzanie prawej i lewej strony alignmentu
                        if len(set(binary_score_right_end)) == 1 and binary_score_right_end[0] == 1: #Sprawdza czy wszystkie elementy listy są takie same za pomocą setu 
                            #(set zwiera tylko unikalne wartości, więc jeśli elementy listy były takie same powinien być tylko 1)
                            
                            #Pętla od początku
                            for k in indexer:
                                if k != i:
                                    seq3 = str(cons_seqs_df.at[k, 'Sequence'])
                                    alignment2 = str(pairwise2.align.globalxx(seq1, seq3, one_alignment_only=True))
                                    search4 = re.search('(?<=seqB=\')[ACTGX\-]+\'', alignment)
                                    search5 = re.search('(?<=score=)[0-9]+', alignment)
                                    search6 = re.search('(?<=end=)[0-9]+', alignment)
                                    aligned_seq_l = search4.group(0)
                                    score_l = int(search5.group(0))
                                    align_length_l = int(search6.group(0))
                                
                                    identity_l = score/align_length_l*100
                                    if identity_l < identity_threshold:
                                        seq2_list_l = list(aligned_seq_l)
                                        # Tworzy binary score
                                        binary_score_l = []
                                        for base in seq2_list_l:
                                            if base == 'A' or base == 'T' or base == 'C' or base == 'G':
                                                binary_score_l.append(1)
                                            else:
                                                binary_score_l.append(0)
                                            
                                        binary_score_left_end2 = binary_score_l[-min_chimera_length:]
                                    
                                        if len(set(binary_score_left_end2)) == 1 and binary_score_left_end2[0] == 1:
                                            
                                            #Sekwencja jest chimerą
                                            reversed_df.at[j, 'Chimera'] = 'chimera'
                            
                            
                        elif len(set(binary_score_left_end)) == 1 and binary_score_left_end[0] == 1:
                            #Pętla od początku
                            #Pętla od początku
                            for k in indexer:
                                if k != i:
                                    seq3 = str(cons_seqs_df.at[k, 'Sequence'])
                                    alignment2 = str(pairwise2.align.globalxx(seq1, seq3, one_alignment_only=True))
                                    search4 = re.search('(?<=seqB=\')[ACTGX\-]+\'', alignment)
                                    search5 = re.search('(?<=score=)[0-9]+', alignment)
                                    search6 = re.search('(?<=end=)[0-9]+', alignment)
                                    aligned_seq_r = search4.group(0)
                                    score_r = int(search5.group(0))
                                    align_length_r = int(search6.group(0))
                                
                                    identity_r = score/align_length_r*100
                                    if identity_r < identity_threshold:
                                        seq2_list_r = list(aligned_seq_r)
                                        # Tworzy binary score
                                        binary_score_r = []
                                        for base in seq2_list_r:
                                            if base == 'A' or base == 'T' or base == 'C' or base == 'G':
                                                binary_score_r.append(1)
                                            else:
                                                binary_score_r.append(0)
                                            
                                        binary_score_right_end2 = binary_score_r[-min_chimera_length:]
                                    
                                        if len(set(binary_score_right_end2)) == 1 and binary_score_right_end2[0] == 1:
                                            
                                            #Sekwencja jest chimerą
                                            reversed_df.at[j, 'Chimera'] = 'chimera'
            
            
            # Filtrowanie po chimerach, przywracanie właściwej kolejności po głębi i przywrócenie właściwych indeksów
            filtered_seqs = reversed_df.drop(reversed_df[reversed_df.Chimera == 'chimera'].index)
            filtered_seqs.sort_values(['Depth'], inplace=True)
            cons_seqs_df = filtered_seqs.reset_index(drop=True)
            
#################### END OF CLASS FUNCTIONS ##########################
######################################################################  
    
################### FUNCTIONS #################################
###############################################################

def new_method():
    
    prework = Prework()
    
    #Ekstrakcja sekwencji z gz
    for filename in os.scandir(pre_input):
        os.system('gunzip ' + filename.path)
    
    # Utworzenie listy z nazwami plików
    all_files_name = os.listdir(pre_input)
    all_files_name.sort()
    
    left_files = all_files_name[::2]
    right_files = all_files_name[1::2]
    
    names_list = []
    pattern = r'_R1_001'
    for name in left_files:
        name_stripped = re.sub(pattern, '', name)
        names_list.append(name_stripped)
    
    # Ilość osobników analizowanych
    iterator = [i for i in range(len(left_files))]
    
    # Merge plików
    for i in iterator:
        left_file_name = left_files[i]
        left_file = f'{pre_input}/{left_file_name}'
        right_file_name = right_files[i]
        right_file = f'{pre_input}/{right_file_name}'
        file_name = names_list[i]
        
        prework.pear(left_file, right_file, file_name, merged, threads)
        
    current = os.getcwd()
    os.chdir(merged)
    
    os.system('ls | grep -v assembled.fastq | xargs rm')
    
    os.chdir(current)
    
    # Oczyszczenie plików
    for filename in os.scandir(merged):
        name = os.path.basename(filename.path)
        prework.grep(filename.path, name, trimmed)
        
    # Usuwanie starterów
    primers = pd.read_csv(primers_file)
    primer_f = primers.at[0, 'primer_f']
    primer_r_rc = str(Seq(primers.at[0, 'primer_r']).reverse_complement())
    region_f = len(primer_f)+1
    region_r = (len(primer_r_rc)+1)*-1
    
    for filename in os.scandir(trimmed):
        name = os.path.basename(filename.path)
        prework.cut_primers(filename.path, name, primer_f, primer_r_rc, region_f, region_r, clean)
        
        
        
def old_method():
    
    prework = Prework()
    
    #Ekstrakcja sekwencji z gz
    for filename in os.scandir(pre_input):
        os.system('gunzip ' + filename.path)
      
    #Utworzenie listy z nazwami plików
    all_files_name = os.listdir(pre_input)
    all_files_name.sort()
    
    left_file_name = all_files_name[0]
    right_file_name = all_files_name[1]
    
    pattern = r'_R1_001'
    name_stripped = re.sub(pattern, '', left_file_name)

    left_file = f'{pre_input}/{left_file_name}'
    right_file = f'{pre_input}/{right_file_name}'
    
    prework.pear(left_file, right_file, name_stripped, merged, threads)
    
    current = os.getcwd()
    os.chdir(merged)
    
    os.system('ls | grep -v assembled.fastq | xargs rm')
    
    os.chdir(current)
    
    # Oczyszczanie plików
    for filename in os.scandir(merged):
        name = os.path.basename(filename.path)
        prework.grep(filename.path, name, trimmed)
        
    ### Usuwanie starterów
    
    # Wczytanie pliku barkodów
    data = pd.read_csv(barcode_file)
    
    # Wyznaczenie starterów
    primer_f = data.at[0, 'primer_f']
    primer_r_rc = str(Seq(data.at[0, 'primer_r']).reverse_complement())
    region_f = len(primer_f)+1
    region_r = (len(primer_r_rc)+1)*-1
    
    iterator = [i for i in range(2, len(data))]
    
    name_file = None
    for filename in os.scandir(trimmed):
        name_file = filename.path
        
    records = list(SeqIO.parse(name_file, "fastq"))
    
    all_amplicons = pd.DataFrame(columns = ['Amplicon', 'ID', 'Sequence'])
    
    len_rec = len(records)
    indexer = [i for i in range(len_rec)]
    
    give = r'{0,3}'
    
    for j in indexer:
        checked_seq = records[j].seq
        checked_id = str(records[j].id)
        reversedcom_seq = checked_seq.reverse_complement()
        for i in iterator:
            amplicon_df = pd.DataFrame(columns = ['Amplicon', 'ID', 'Sequence'])
            barA = data.at[i, 'primer_f']
            barB = data.at[i, 'primer_r']
            patternA = f'\A[A-Z]{give}{barA}'
            patternB = f'{barB}[A-Z]{give}$'
            if re.search(patternA, str(checked_seq)) != None and re.search(patternB, str(checked_seq)) != None:
                strippedA = re.sub(patternA, '', str(checked_seq))
                strippedB = re.sub(patternB, '', strippedA)
                amplicon_df.at[0, 'Amplicon'] = data.at[i, 'sample']
                amplicon_df.at[0, 'ID'] = checked_id
                amplicon_df.at[0, 'Sequence'] = strippedB
            elif re.search(patternA, str(reversedcom_seq)) != None and re.search(patternB, str(reversedcom_seq)) != None:
                strippedA = re.sub(patternA, '', str(reversedcom_seq))
                strippedB = re.sub(patternB, '', strippedA)
                amplicon_df.at[0, 'Amplicon'] = data.at[i, 'sample']
                amplicon_df.at[0, 'ID'] = checked_id
                amplicon_df.at[0, 'Sequence'] = strippedB
                
            #all_amplicons = all_amplicons.append(amplicon_df.iloc[0])
            frames = [all_amplicons, amplicon_df]
            all_amplicons = pd.concat(frames, ignore_index=True)
            
    all_amplicons.sort_values(by=['Amplicon'], inplace=True)
    
    # Tworzenie fasta dla poszczególnych amplikonów
    for k in iterator:
        amplicon_number = data.at[k, 'sample']
        amplicon_f_df = all_amplicons[all_amplicons['Amplicon'] == amplicon_number]
        outpath = f'{disbarred_dir}/amplicon_{amplicon_number}.fa'
        with open(outpath, 'w') as outfile:
            for index, row in amplicon_f_df [ ['ID', 'Sequence'] ].iterrows():
                outfile.write('>' + row.ID + '\n' + row.Sequence + '\n')
                
        # Wycięcie starterów
        command = f'cat {outpath} | seqkit amplicon -F {primer_f} -R {primer_r_rc} -r {region_f}:{region_r} -o {clean}/amplicon_{amplicon_number}.fa'
        process = subprocess.check_output(command, shell=True)



def fazeI(pair):
    
    files, filepath = pair
    
    go = Amplicon(files)
    go.parse_sequence_file(filepath, method_end)
        
    #Odrzucenie amplikonów o zbyt małej głębokości
    depth_of_amplicon = go.amplicon_table['Depth'].sum()
    if depth_of_amplicon >= min_amplicon_depth:
        
        Stats_table_amplicon_depth = pd.DataFrame(columns = ['Amplicon', 'Amplicon_updated', 'Depth'])
        Stats_table_amplicon_depth.at[files-1, 'Depth'] = depth_of_amplicon
        Stats_table_amplicon_depth.at[files-1, 'Amplicon'] = files
            
        go.find_minmax()
        go.create_fasta(p7)
        
        seqs_fasta_dir = f'{p8}/amplicon_{files}'
        os.mkdir(seqs_fasta_dir)
        go.create_seq_fasta(seqs_fasta_dir) #Tworzy fasty dla każdej sekwencji w amplikonie
        
        clusters_dir = f'{p9}/amplicon_{files}'
        os.mkdir(clusters_dir)
        
        shutil.copy('./fogsaa',clusters_dir)
        current = os.getcwd()
        os.chdir(clusters_dir)
        
        do = Cluster(files)   #Wczytanie pliku z funkcjami klastrowania
        
        all_clusters = pd.DataFrame(columns = ['ID', 'Sequence', 'Depth', 'Length'])
        
        go.amplicon_table['state'] = 'U' #Dodanie kolumny 'state'
        
        index_list = list(go.amplicon_table.index.values)
        
        for index in index_list:
            if go.amplicon_table.at[index, 'state'] == 'U':
                if go.amplicon_table.at[index, 'Frequency'] >= min_amplicon_seq_frequency:
                    cluster_pd = pd.DataFrame(columns = ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'])
            
                    dom_seq_row = go.amplicon_table.iloc[index]
                    dom_seq_hash = go.amplicon_table.at[index, 'Hash']
                    dom_seq_length = go.amplicon_table.at[index, 'Length']
                    dom_seq_depth = go.amplicon_table.at[index, 'Depth']
            
                    seq1 = f'../../amplicons_filtered_seqs/amplicon_{files}/{dom_seq_hash}.fasta'
            
                    cluster_pd = cluster_pd.append(dom_seq_row)
                    
                    go.amplicon_table.at[index, 'state'] = 'C'
                    
                    index_klastra = index
            
                    for index in index_list:
                        if go.amplicon_table.at[index, 'state'] == 'U':
                            if go.amplicon_table.at[index, 'Length'] == dom_seq_length:
                                seq2_hash = go.amplicon_table.at[index, 'Hash']
                                seq2 = f'../../amplicons_filtered_seqs/amplicon_{files}/{seq2_hash}.fasta'          
                                do.seq2seq(seq1, seq2, substitution_error_threshold, index)
                                                    
                                if substitution_error < substitution_error_threshold:
                                    freq_to_dom = go.amplicon_table.at[index, 'Depth']/dom_seq_depth*100
                        
                                    if freq_to_dom < min_dominant_frequency_threshold:

                                        cluster_pd = cluster_pd.append(go.amplicon_table.loc[index]) #dodaje sekwencję do klastra dominującego
                                        go.amplicon_table.at[index, 'state'] = 'C'

                    do.multiple_aline_seqs(cluster_pd, "./", index_klastra)
                               

                    all_clusters = all_clusters.append(consensus_df.iloc[0])
        
        #Stworzenie ostatecznego pliku fasta z sekwencjami konsensusowymi dla danego amplikonu
        
        outpath = './consensus_seqs.fasta'
        with open(outpath, 'w') as outfile:
            for index, row in all_clusters[ ['ID', 'Sequence', 'Depth', 'Length'] ].iterrows():
                outfile.write('>' + row.ID + ' | depth: ' + str(row. Depth) + ' | length: ' + str(row.Length) + '\n' + row.Sequence + '\n')
                
        # Stworzenie pliku csv z sewkencjami konsenusowymi
        outpath2 = './consensus_seqs.csv'
        all_clusters.to_csv(outpath2, index=False)
        os.chdir(current)

        
        outpath3 = f'{clusters_dir}/amplicon_info.csv'
        Stats_table_amplicon_depth.to_csv(outpath3, index=False)
                        
################ END OF FUNCTIONS ########################
##########################################################

#Prework start
if method == 'new':
	method_end = 'fastq'
	new_method()
elif method == 'old':
	method_end = 'fasta'
	old_method()

	


### Stworzenie listy z numerami i ścieżkami do plików dla klastrowania

files = 0
command_list = []
for filename in os.scandir(clean):
    if filename.is_file():
        files = files+1
        path = str(filename.path)
        line = files, f'{path}'
        command_list.append(line)


       
def fazeII():
    start = time.time()
    pool = Pool(processes=4)
    result = pool.map(fazeI, command_list)
    pool.close()
    end = time.time()
    delta = end-start
    print(f'Operacja zabrała {delta:.3f} sekund')
    
fazeII()
        

### FILTROWANIE ###


all_allels = pd.DataFrame(columns = ['Amplicon', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'])

# Ilość amplikonów
files = 0
for filename in os.scandir(p9):
    if filename.is_dir():
        files = files+1

# Lista ilości amplikonów        
no_files = [i for i in range(1, files+1)]


for i in no_files:
    con_seqs_dir = f'{p9}/amplicon_{i}/consensus_seqs.csv'
    cons_seqs_df = pd.read_csv(con_seqs_dir)
    be = Filter()
    
    depth_of_amplicon = cons_seqs_df['Depth'].sum()  # Wylicza głębię całego amplikonu dla sekwencji konsensusowych
    cons_seqs_df['Frequency'] = cons_seqs_df['Depth']/depth_of_amplicon*100 # Wylicza częstość występowania sekwencji w amplikonie
    
    #Odrzuca sewkencje o zbyt niskiej głębi
    cons_seqs_df.drop(cons_seqs_df[cons_seqs_df.Depth < min_allel_depth].index, inplace=True)
    #Odrzuca sekwencje o zbyt niskiej czestotliwości
    cons_seqs_df.drop(cons_seqs_df[cons_seqs_df.Frequency < min_amplicon_seq_frequency].index, inplace=True)
    
    cons_seqs_df.reset_index(drop=True, inplace=True)

    
    # CHIMERY
    be.is_chimera(cons_seqs_df)

    
    #Odrzucanie nadmiernej ilości alleli w amplikonie
    seqs_count = len(cons_seqs_df)
    if seqs_count > max_allele_number:
        cons_seqs_df = cons_seqs_df[:max_allele_number]

    
    #Dodanie kolumny z nazwą amplikonu
    cons_seqs_df.insert(0, 'Amplicon', i)
    
    #Dodanie linijki
    all_allels = all_allels.append(cons_seqs_df)
    
    #Głębia amplikonu
    a_info_dir = f'{p9}/amplicon_{i}/amplicon_info.csv'
    a_info_df = pd.read_csv(a_info_dir)
    current_depth = a_info_df.at[0, 'Depth']
    
    #Stworzenie ostatecznego pliku fasta z allelami dla danego amplikonu
        
    outpath = f'{p10}/amplicon{i}.fasta'
    with open(outpath, 'w') as outfile:
        outfile.write('Amplicon depth: ' + str(current_depth) + '\n')
        for index, row in cons_seqs_df[ ['ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
            outfile.write('>' + row.ID + ' | depth: ' + str(row. Depth) + ' | length: ' + str(row.Length) + ' | frequency: ' + str(row.Frequency) + '\n' + row.Sequence + '\n')
    
all_allels.reset_index(drop=True, inplace=True)
whole_depth =  all_allels['Depth'].sum()  

droped_duplcates = all_allels.drop_duplicates(subset=['Sequence'])
labels = droped_duplcates['Sequence'].tolist()

#all_allels.groupby(labels)

final_allels = pd.DataFrame(columns = ['Amplicon', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency', 'Frequency_across_amplicons'])
final_stats = pd.DataFrame(columns = ['Sequence', 'Amplicons', 'Count', 'Depth', 'Frequency_across_amplicons', 'Mean_freq', 'Max_freq', 'Min_freq'])

# len_df = len(all_allels)
# iterator = [i for i in range(len_df)]

for label in labels:
    allel_df = pd.DataFrame(columns = ['Amplicon', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency', 'Frequency_across_amplicons'])
    allel_stats = pd.DataFrame(columns = ['Sequence', 'Amplicons', 'Count', 'Depth', 'Frequency_across_amplicons', 'Mean_freq', 'Max_freq', 'Min_freq'])
    
    for index, row in all_allels[ ['Amplicon', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
#     for x in iterator:
        if all_allels.at[index, 'Sequence'] == label:
        #if row.Sequence == label:
            allel_df = allel_df.append(all_allels.iloc[index])
    
    allel_df.reset_index(drop=True, inplace=True)
    allel_depth = allel_df['Depth'].sum() #Głębia allela na przestrzeni wszystkich amplikonów
    allel_df['Frequency_across_amplicons'] = allel_depth/whole_depth*100 #Częstość występowania na przestrzeni wszystkich ampikonów
    count = len(allel_df) #Ilość wystąpień danego allela we wszustkich amplikonach
    
    amplicon_list = [] 
    for index, row in allel_df[ ['Amplicon', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency', 'Frequency_across_amplicons'] ].iterrows():
        amplicon_list.append(row.Amplicon)
    
    allel_stats.at[0, 'Sequence'] = allel_df.at[0, 'Sequence']
    allel_stats.at[0, 'Amplicons'] = str(amplicon_list)
    allel_stats.at[0, 'Count'] = count
    allel_stats.at[0, 'Depth'] = allel_depth
    allel_stats.at[0, 'Frequency_across_amplicons'] = round(allel_depth/whole_depth*100, 2)
    allel_stats.at[0, 'Mean_freq'] = round(allel_df['Frequency'].mean(), 2)
    allel_stats.at[0, 'Max_freq'] = round(allel_df['Frequency'].max(), 2)
    allel_stats.at[0, 'Min_freq'] = round(allel_df['Frequency'].min(), 2)
    
    final_stats = final_stats.append(allel_stats)
    final_stats.sort_values(['Frequency_across_amplicons'], ascending=False, inplace=True)
    
    #Dołącznie allela z częstością do df końcowego
    final_allels = final_allels.append(allel_df)
    final_allels.sort_values(['Amplicon'], inplace=True)
    
# Stworzenie outputu fasta dla ostatecznych sekwencji
outpath = f'{p10}/all_allels.fasta'
with open(outpath, 'w') as outfile:
    for index, row in final_allels[ ['Amplicon', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency', 'Frequency_across_amplicons'] ].iterrows():
        outfile.write('>' + 'Amplicon: ' + str(row.Amplicon) + ' | ' + row.ID + ' | depth: ' + str(row. Depth) + ' | length: ' + str(row.Length) + ' | frequency: ' + str(row.Frequency) + ' | frequency across amplicons: ' + str(row.Frequency_across_amplicons) + '\n' + row.Sequence + '\n')

#Stworzenie outputu tsv dla ostatecznych sekwencji
outpath = f'{p10}/all_allels_stats.csv'
final_stats.to_csv(outpath, index=False)


###########################################################################################################################################
####################################################### THAT'S ALL FOLKS ##################################################################











