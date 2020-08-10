
# coding: utf-8

# # Gene Annotation Project
# In this project, we are primarily concerned with annotating the ebony gene from flyBase. The end goal of this python program is to anootate the ebony gene in the SnapGene viewer. 

# In[1]:


from Bio import SeqIO
import pandas as pd
from datetime import datetime as dt

#extract sequence from FASTA file
def openFASTA(file_name):
    for seq_record in SeqIO.parse(file_name + ".fasta", "fasta"):
        return seq_record.seq
        
def start_stop(seq_1, seq_2):
    if len(seq_1) > len(seq_2):
        main = seq_1
        target = seq_2
    else: 
        main = seq_2
        target = seq_1
    if main.count(target) == 1:
        loc_1 = main.find(target) + 1
        loc_1_final = loc_1 - 1
        loc_2 = loc_1_final + len(target)
        start_end = [loc_1, loc_2]
        return start_end 
    else: 
        return 'Error: Multiple instances of this exon detected'  
    
def open_multi_FASTA(file_name):
    from Bio import SeqIO
    exons = {}
    for seq_record in SeqIO.parse(file_name + ".fasta", "fasta"):
        gene_id = seq_record.id
        exons[gene_id] = seq_record.seq
    return exons

open_multi_FASTA('all_files')


# In[2]:


main_file = input ('Enter name of the extended gene region file: ')
target_file =  input("Enter name of exon file: ") 
sequence_type = input('Enter Sequence Type: ')
multiple_files = input('Are the region(s) represented in this file continuous? (Input Y or N): ' )

first_terms = []
last_terms = []
final_dict = {}
exon_names = []
cols = [sequence_type, 'Start Nucleotide Position', 'Stop Nucleotide Position']

ext_region = openFASTA(main_file)

if multiple_files == 'N':
    exons = open_multi_FASTA(target_file)
    for gene_id in exons:
        seq = exons[gene_id]
        point = start_stop(seq, ext_region)
        final_dict[gene_id] = point

else:
    exons = openFASTA(target_file)
    seq = exons[gene_id]
    point = start_stop(seq, ext_region)
    final_dict[gene_id] = point


start_stop_df = pd.DataFrame(final_dict)
start_stop_df = start_stop_df.transpose()
start_stop_df.reset_index(inplace = True, level = 0)
start_stop_df.columns = ['exon_id','start_nucl_pos', 'stop_nucl_pos'] 
start_stop_df = start_stop_df.sort_values(by = ['start_nucl_pos'], ascending = True)

today = dt.today().strftime('%Y-%m-%d')
output_file = 'ebony_gene_exons_' + today + '.csv'
start_stop_df.to_csv(output_file, index = False)

display(start_stop_df)

