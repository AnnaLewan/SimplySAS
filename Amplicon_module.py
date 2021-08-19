#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import pandas as pd
import uuid
#import re
import subprocess

class Amplicon:
    def __init__(self, number):
        self.number = number
        self.clusters = None
        self.seqs_count = None
        
        
    def cut_primers(self, input_fastq, outputs_dir_good, outputs_dir_bad, starter_f, starter_r):
        command = f'python3 cutPrimers.py -r1 {input_fastq} -pr15 {starter_f} -pr13 {starter_r} -tr1 {outputs_dir_good}/amplicon{self.number}_primerfree.fastq -utr1 {outputs_dir_bad}/odrzucone_ampli{self.number}.fastq'        
        process = subprocess.check_output(command, shell=True)
        
            
    def parse_sequence_file(self, inputfastq_primerfree):
        global amplicon_table
        amplicon_table = pd.DataFrame(columns = ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'])      ## ogarnąc czy przypadkiem nie trzeba zrobić w inicie  
        
        dedup_records = defaultdict(list) #tworzy słownik, dla których wartościami jest lista

        for record in SeqIO.parse(inputfastq_primerfree, "fastq"):
            dedup_records[str(record.seq)].append(record.id) #tworzy słownik, dla którego key=seq, a value=id. Na tym etapie dochodzi do odnalezienia duplikatów - jest ich tyle wartości id
        for seq, ids in sorted(dedup_records.items(), key=lambda t: len(t[1]), reverse=True): #ta linijka służy przesortowaniu sekwencji pod względem ich głębokości                
            id_hold = ids[0]
            id_len = len(ids)
            seq_len = len(seq)
            unique_hash = str(uuid.uuid1())
            amplicon_table = amplicon_table.append({'Hash' : unique_hash, 'ID' : id_hold, 'Sequence' : seq, 'Depth' : id_len, 'Length' : seq_len}, ignore_index = True)
        
        self.seqs_count = len(amplicon_table) # Ilość sekwencji w amplikonie
        depth_of_amplicon = amplicon_table['Depth'].sum() #Wylicza głębie całego amplikonu
        
        for index, row in amplicon_table[ ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
            freq_in_amplicon = row.Depth/depth_of_amplicon*100
            row.Frequency_per_amplicon = freq_in_amplicon            
        
        return amplicon_table

    
    def create_fasta(self, inputs_dir):
        outpath = f'{inputs_dir}/{self.number}_amplicon.fasta'
        with open(outpath, 'w') as outfile:
            for index, row in amplicon_table[ ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
                outfile.write('>' + row.Hash + ' | ' + row.ID + ' | depth: ' + str(row.Depth) + ' | length: ' + str(row.Length) + ' | frequency per amplicon: ' + str(row.Frequency) + '\n' + row.Sequence + '\n')
                
                
    def create_seq_fasta(self, seqs_fasta_dir):
        for index, row in amplicon_table[ ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
            outpath = f'{seqs_fasta_dir}/{row.Hash}.fasta'
            with open(outpath, 'w') as outfile:
                outfile.write('>' + row.Hash + ',' + row.ID + ',' + str(row.Depth) + ',' + str(row.Length) + ',' + str(row.Frequency) + '\n' + row.Sequence + '\n')


