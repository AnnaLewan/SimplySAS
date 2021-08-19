#!/usr/bin/env python3

import sys
import os
import argparse
import pandas as pd
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align import AlignInfo
import zipfile
import tarfile
import Amplicon_module
import Cluster_module
from collections import defaultdict
import uuid
import subprocess
import shutil

parser = argparse.ArgumentParser(prog='PartSAS', description="Performs part of genotyping of amplicon sequencing data by clustering errors and filtering artefacts following a 3-step based analysis: de-multiplexing, clustering and filtering.")

parser.add_argument("-i",'--INPUT',type=str,required=True,help="Input set of FASTQ or FASTA files packed into a unique .ZIP or .TAR.GZ file.")
parser.add_argument("-d",'--DATA',type=str,required=True,help="CSV file with primer/amplicon data.")
parser.add_argument("-o",'--OUTPUT',type=str,required=True,help="Output folder name.")
parser.add_argument("-ml", '--MINLEN',type=int,required=True,help="Minimal length of sequence to consider clustering.")
parser.add_argument("-el", '--EXPLEN',type=int,required=True,help="Expected length of the marker sequence.")
parser.add_argument("-se", '--SUBERROR',type=float,required=True,help="Threshold for permissible substitution error (%).")
parser.add_argument("-df", '--DOMFREQ',type=float,required=True,help="Minimum frequency respect to the dominant (%).")
parser.add_argument("-af", '--AMFREQ',type=float,required=True,help="Minimum frequency per-amplicon (%).")

args = parser.parse_args()

input_multifile = args.INPUT
input_csv = args.DATA
output_dir = args.OUTPUT
min_length = args.MINLEN #minimalna długość sekwencji aby wziąć ją pod uwagę w czasie klastrowania
expected_len = args.EXPLEN #oczekiwana długość markera
substitution_error_threshold = args.SUBERROR #threshold dla dopuszczalnego błędu substytucji aby odczyt mógł zostać przypisany do danego klastra (w %)
min_dominant_frequency_threshold = args.DOMFREQ
min_amplicon_seq_frequency = args.AMFREQ


p1 = f'{output_dir}/ampliSAS_analysis'
os.makedirs(p1)

p2 = f'{output_dir}/ampliSAS_analysis/input'
os.makedirs(p2)

if tarfile.is_tarfile(input_multifile) or zipfile.is_zipfile(input_multifile) != True:
    raise Exception('File is not a multifile')
            
if os.path.isfile(input_csv) != True:
    raise Exception('CSV file not existing')
    
primer_file = pd.read_csv(input_csv, names=[0,1])

forward = primer_file.iloc[1, 1]
reverse = str(Seq(primer_file.iloc[2, 1]).reverse_complement())
        
ID = primer_file.iloc[0, 0]
forward_ID = f'{ID}_forward'
reverse_ID = f'{ID}_reverse'

p3 = f'{p1}/input/fasta_primers'
os.mkdir(p3) 
   
outpath_f = f'{p3}/starter_f.fasta'
outpath_r = f'{p3}/starter_r.fasta'
with open(outpath_f, 'w') as outfile:
    outfile.write('>' + forward_ID + '\n' + forward)
with open(outpath_r, 'w') as outfile:
    outfile.write('>' + reverse_ID + '\n' + reverse)
    
print(f'Extracting fastqs from multifile')

#for i in input_multifile:
with zipfile.ZipFile(input_multifile, 'r') as zip:
    zip.extractall(path = f'{p2}')
        
print(f'Extraction finished')

fastqs_folder = re.sub('\.zip', '', input_multifile)

os.rename(f'{p2}/{fastqs_folder}', f'{p2}/fastqs')

p4 = f'{p2}/fastqs'
p5 = f'{p2}/fastqs_primerfree'
p6 = f'{p2}/fastqs_primerfree/out'
os.mkdir(p5)
os.mkdir(p6)

#Wyciananie starterów
files = 0
for filename in os.scandir(p4):
    files = files+1
    if filename.is_file():
        go = Amplicon_module.Amplicon(files)
        go.cut_primers(filename.path, p5, p6, outpath_f, outpath_r)
        
shutil.rmtree(p6)

p7 = f'{p2}/fasta'
os.mkdir(p7)

p8 = f'{p1}/amplicons_filtered_seqs' # Tu będa znajdowały się foldery dla każdego amplikonu zawierające każdą sekwencję w amplikonie jako osobny plik fasta. Sekwencje te są już przefiltrowane względem długości, głebi itd
os.mkdir(p8)

p9 = f'{p1}/clusters'
os.mkdir(p9)

p10 = f'{p1}/clusters/mafftout'
os.mkdir(p10)

files = 0
for filename in os.scandir(p5):
    files = files+1
    if filename.is_file():
        go = Amplicon_module.Amplicon(files)
        go.parse_sequence_file(filename.path)
        go.create_fasta(p7) #SPRAWDZONE
        do = Cluster_module.Cluster(files, min_length, expected_len)   #Wczytanie pliku z funkcjami klastrowania
        do.find_dominant()
        seqs_fasta_dir = f'{p8}/amplicon_{files}'
        os.mkdir(seqs_fasta_dir)
        go.create_seq_fasta(seqs_fasta_dir)
        clusters_dir = f'{p9}/amplicon_{files}'
        os.mkdir(clusters_dir)
        mafftouts = f'{p10}/amplicon_{files}'
        os.mkdir(mafftouts)
        all_clusters = pd.DataFrame(columns = ['ID', 'Sequence', 'Depth', 'Length'])
        subdominants = 0 #zmienna do pilnowania liczby klastrów poza dominującym
        subdominants_df =  pd.DataFrame(columns = ['Subdominant', 'Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency']) #Tworzy dataframe, w którym będzie zapisywana każda sekwencja subdominująca z danego ampilkonu
        
        #Pierwsze przejście - dla sekwencji domminującej
        for index, row in amplicon_table[ ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
            #if index != 0: #po dodaniu dropa w find_dominant() niepotrzebne
                seq1a = cluster_dom.loc[0, 'Hash']
                seq1 = f'{seqs_fasta_dir}/{seq1a}.fasta'
                seq2 = f'{seqs_fasta_dir}/{row.Hash}.fasta'           
                do.seq2seq(seq1, seq2, substitution_error_threshold)
                if substitution_error < substitution_error_threshold:
                    #Warunek freq
                    freq_to_dom = row.Depth/cluster_dom.loc[0, 'Depth']*100
                    if row.Frequency < min_amplicon_seq_frequency or freq_to_dom < min_dominant_frequency_threshold:
                        if row.Length == cluster_dom.loc[0, 'Length']:
                            cluster_dom = cluster_dom.append(amplicon_table.iloc[index]) #dodaje sekwencję do klastra dominującego
                            amplicon_table.drop(index) #usuwa sekwencję z głównego dataframe amplikonu
                            do.multiple_aline_seqs(cluster_dom, clusters_dir, mafftouts, 'dom')
                    else:
                        subdominants += 1
                        subdominants_df = subdominants_df.append(amplicon_table.iloc[index])
                        subdominants_df = subdominants_df.append({'Subdominant' : subdominants}, ignore_index = True)
                        amplicon_table.drop(index)
        all_clusters = all_clusters.append(consensus_df.iloc[0])
        # Na tym kończy się pierwsze przejście przez amplikon. Jego efektem jest wytworzenie klastra,
        # który jako podstawe ma sekwencję dominującą w amplikonie oraz usunięcie z amplicon_table
        # zarówno sekwencji dominującej jak i pozostałych należących do klastra sekwencji.
        # Dodatkowo tworzony jest DataFrame, który przechowuje wszystkie sekwencje, które zostały
        # uznane za subdominujące i mogą zostac podstawami nowych klastrów
        
#         for index, row in subdominants_df[ ['Subdominant', 'Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
#             dom_sub_seq_cluster = pd.DataFrame(columns = ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'])
#             dom_sub_seq = subdominants_df.iloc[0]
#             dom_sub_seq_cluster = dom_sub_seq_cluster.append( [dom_sub_seq] )
#             seq1a = dom_sub_seq.loc[0, 'Hash']
#             seq1 = f'{seqs_fasta_dir}/{seq1a}.fasta'
#             for index, row in subdominants_df[ ['Subdominant', 'Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
#                 if index != 0:
#                     if row.Length == dom_sub_seq.loc[0, 'Length']:
#                         seq2 = f'{seqs_fasta_dir}/{row.Hash}.fasta'
#                         do.seq2seq(seq1, seq2, substitution_error_threshold)
#                         if substitution_error < substitution_error_threshold:
#                             dom_sub_seq_cluster = dom_sub_seq_cluster.append(subdominants_df.iloc[index])
#                             subdominants_df.drop(subdominants_df.iloc[index])
#                             subdominants_df.reset_index(drop=True, inplace=True)
#             for index, row in amplicon_table[ ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
#                 dom_sub_seq = dom_sub_seq_cluster.iloc[0]
        
        #Stworzenie ostatecznego pliku fasta z sekwencjami konsensusowymi dla danego klastra
        
        outpath = f'{clusters_dir}/consensus_seqs.fasta'
        with open(outpath, 'w') as outfile:
            for index, row in all_clusters[ ['ID', 'Sequence', 'Depth', 'Length'] ].iterrows():
                outfile.write('>' + row.ID + ' | depth: ' + str(row. Depth) + ' | length: ' + str(row.Length) + '\n' + row.Sequence + '\n')
        



