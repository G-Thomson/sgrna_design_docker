#! /usr/bin/env python

import os
import subprocess
import shutil

import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from getGene import getGene

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument('geneList',
                    type=str,
                    help='')
args = parser.parse_args()

# Load data
df = pd.read_csv(args.geneList)

# Make exon files
file_path_parent = '/home/Pipeline/'
os.mkdir(file_path_parent+'tmp_indiv_files/')

df_grouped = df.groupby('orthogroup')

#Loop over each orthogroup
for group_name, df_group in df_grouped:
    # group_name - orthogroup
    # df_group - dataframe
    file_path = file_path_parent+'tmp_indiv_files/'+group_name+'.fa'

    with open(file_path, mode='w') as file_object:

        # Loop over each gene
        for grp_ind, grp_df in df_group.iterrows():
            transcript = grp_df[0]
            cds = getGene(transcript, mode = 'exons', genome = '/home/Data/C.sinensis_v2.0_HZAU_csi.chromosome.fa', database = "/home/Data/citrus.db")

            SeqIO.write(cds, file_object, "fasta")
        
    file_object.close() 


import itertools
for group_name, df_group in itertools.islice(df_grouped, 177, 180): #<<<<<<<<<
# for group_name, df_group in df_grouped: 
    # # group_name - orthogroup
    # # df_group - dataframe
    file_path = file_path_parent+'tmp_indiv_files/'+group_name+'.fa'
    os.mkdir(file_path_parent+'crispys_output/')

    crispys = "python Stage0.py "+file_path+" "+file_path_parent+'crispys_output/'
    subprocess.call(crispys, shell=True)

    print("CHECKPOINT 1")

    if not os.path.exists(file_path_parent+'crispys_output/output.csv'):
        shutil.rmtree(file_path_parent+'crispys_output/')
        continue

    print("CHECKPOINT 2")

    R_cleanup = "Rscript "+file_path_parent+"/Format_sgRNA_table.R --in "+file_path_parent+'crispys_output/output.csv'
    subprocess.call(R_cleanup, shell=True)

    print("CHECKPOINT 3")

    if not os.path.exists("Formated_sgRNA_table.csv"):
        shutil.rmtree(file_path_parent+'crispys_output/')
        continue

    print("CHECKPOINT 4")

    sgRNA_call="python "+file_path_parent+"sgRNA_selection.py Formated_sgRNA_table.csv "+group_name
    subprocess.call(sgRNA_call, shell=True)

    print("CHECKPOINT 5")

    os.remove("Formated_sgRNA_table.csv")
    shutil.rmtree(file_path_parent+'crispys_output/')

shutil.rmtree(file_path_parent+'tmp_indiv_files/')

