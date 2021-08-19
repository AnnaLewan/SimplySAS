#!/usr/bin/env python3

import pandas as pd
import subprocess
import re
from Bio import AlignIO
from Bio.Align import AlignInfo

class Cluster:
    def __init__(self, amplicon, min_length, expected_len):
        self.dominant_seq = None
        self.amplicon = amplicon
        self.min_length = min_length
        self.expected_len = expected_len
        
        
    def find_dominant(self):
#        amplicon_table.sort_values('Lenght', ascending=False, inplace=True)
        amplicon_table.drop(amplicon_table[amplicon_table.Length < self.min_length].index, inplace=True)  #odrzuca sekwencje o długości niższej niż zadana przez użytkownika
        amplicon_table.reset_index(drop=True, inplace=True)
        
        
        global cluster_dom
        cluster_dom = pd.DataFrame(columns = ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'])
#        if amplicon_table[amplicon_table.Length == self.expected_length].iloc[0]:
#            dominant_seq = amplicon_table.iloc[0]
#        else:
        for index, row in amplicon_table[ ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
            if row.Length == self.expected_len:
                dominant_seq = amplicon_table.iloc[0]
                cluster_dom = cluster_dom.append( [dominant_seq] )
                amplicon_table.drop(0)
                amplicon_table.reset_index(drop=True, inplace=True)
                break
                       
                
    def seq2seq(self, seq1, seq2, substitution_error_threshold):
        command = f'./fogsaa {seq1} {seq2} 1 0 1 -1 -1'
        process = subprocess.check_output(command, shell=True) #Przykładowy wynik: b'176\n176\nElapsed time: 1 milliseconds\ntotal nodes expanded==176\n\nscore= 174\n'
        process = process.decode('utf-8') #Zamiana bajtów na str
                
        score = int(re.sub('score=\ ', '', ' '.join(map(str, re.findall('score=\ .[0-9]*', process))))) # Wyodrębnienie score dla alignmentu
        len_seqs = int(re.sub('score=\ ', '', ' '.join(map(str, re.findall('(?<=[0-9]\s)[0-9]+', process))))) # Wyodrębnia długość porównywanych sekwencji - założenie, że porównywane sekwencje są tylko tej samej długości
        
        # Duży problem a propos wyliczenia substitution_error. Fogsaa nie zwraca jaka część score to kara za mismatch.
        # Oznacza to, że nie wiem ile substytucji tak naprawdę zaszło, a bez tego nie mogę wyliczyć błędu substytucji.
        
        global substitution_error
        substitution_error = (len_seqs - score)/len_seqs*100 # Wylicza poziom błedów dla danego alignmentu
        
        
    def multiple_aline_seqs(self, cluster_df, clusters_dir, mafftout_dir, number):
        #Tworzenie wspólnej fasty dla klastra
        outpath = f'{clusters_dir}/{number}_cluster.fasta'
        with open(outpath, 'w') as outfile:
            for index, row in cluster_df[ ['Hash', 'ID', 'Sequence', 'Depth', 'Length', 'Frequency'] ].iterrows():
                outfile.write('>' + row.Hash + ' | ' + row.ID + ' | depth: ' + str(row.Depth) + ' | length: ' + str(row.Length) + ' | frequency per amplicon: ' + str(row.Frequency) + '\n' + row.Sequence + '\n')
                        
        #Mafft na pliku --> nie jest potrzebny jesli do klastra trafiają tylko sekwencje o tej samej długości
       # command = f'mafft --auto --quiet {outpath} > {mafftout_dir}'
       # process = subprocess.check_output(command, shell=True)
        
        # Stworzenie dataframe na sekwencję konsensusową z danymi
        global consensus_df
        
        consensus_df = pd.DataFrame(columns = ['ID', 'Sequence', 'Depth', 'Length'])
        id_cluster = f'cluster_{number}'
        len_cluster = cluster_df.loc[0, 'Length']
        depth_cluster = cluster_df['Depth'].sum()
               
        #Ustalenie konsensusowej
        align = AlignIO.read(outpath, 'fasta')
        summary_align = AlignInfo.SummaryInfo(align)
        consensus = str(summary_align.dumb_consensus())
    #    my_pssm = summary_align.pos_specific_score_matrix(consensus,chars_to_ignore = ['N']) # Matrix z częstością dla danej pozycji 
    
        consensus_df = consensus_df.append({'ID' : id_cluster, 'Depth' : depth_cluster, 'Length' : len_cluster, 'Sequence' : consensus}, ignore_index = True)
