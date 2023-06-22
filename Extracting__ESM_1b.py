#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from pandas import DataFrame
from sklearn.preprocessing import StandardScaler, normalize
data = pd.read_excel('../data/5----all_MMdata_pro_seq_simple_label_var_22.xlsx', sheet_name='set_for_ESM1b')
data.head(5)
Sequence = data['Sequence'].tolist()
ENSP = data['ENSP'].tolist()
Name = data['Name'].tolist()
Variant_AA_Change_original = data['Variant_AA_Change_original'].tolist()
len_seq = []
for i in range(len(Sequence)):
    if(Sequence[i]=='-'):
        len_seq.append(0)
    else:
        len_seq.append(len(Sequence[i]))
len(len_seq)
SEQ = Sequence
print(len(ENSP),len(Sequence))

def ensp_name(cut_num,ensp_id):
    lst = []
    for i in range(cut_num):
        tp_name = ensp_id + '_' + str(i)
        lst.append(tp_name)
    return(lst)
# lst = ensp_name(33, 'ENSP00000467141')
# print(lst)
ENSP_keep = []
SEQ_keep = []
for i in range(len(ENSP)):
    temp_seq = SEQ[i]
    temp_ensp = ENSP[i]
    if(len(Sequence[i])>1022):
        ## ENSP
        cut_num = int(len(Sequence[i])/1022 + 1 )
        lst = ensp_name(cut_num,temp_ensp)
        for k in range(len(lst)):
            ENSP_keep.append(lst[k])
        ## SEQ
        for i in range(0, len(Sequence[i]), 1022):
            tp = temp_seq[i:i+1022]
            SEQ_keep.append(tp)
    else:
        ENSP_keep.append(temp_ensp)
        SEQ_keep.append(temp_seq)            
print(len(ENSP_keep),len(SEQ_keep)) 
file = open('../data/5----all_MMdata_4ESM1b_seqs_cutinto_1022AA.txt','w')
for y in range(len(ENSP_keep)):
    file.write('>'+ENSP_keep[y]+'\n')
    file.write(SEQ_keep[y]+'\n')
file.close()
print('save sucess !')

## 2. 提取已经处理过的ENSP名称
import os
import pandas as pd
path ="../data/MissenseM_ALL_ESM1bEmbedder/"
datanames = os.listdir(path)
print(datanames, type(datanames))
processed_name = []
for p in range(len(datanames)):
    temp = datanames[p].split('_')[0]
    processed_name.append(temp)
print(processed_name)
from bio_embeddings.embed import SeqVecEmbedder, ProtTransBertBFDEmbedder, ESM1bEmbedder,ESM1vEmbedder,ProtTransT5XLU50Embedder
from Bio import SeqIO
import pandas as pd
import numpy as np
sequences = []
sequences_names = []
for record in SeqIO.parse("../data/5----all_MMdata_4ESM1b_seqs_cutinto_1022AA.fasta", "fasta"):
    if (record.id not in processed_name):
        sequences.append(record)
        sequences_names.append(record.id)
    else:
        continue
    print('name--seq length:\t', record.id, len(record.seq))

print('length：', len(sequences),len(sequences_names))
embedder = ESM1bEmbedder()
for i in range(0,len(sequences_names),2):
    j = i
    temp_seq = sequences[j:j+2]
    k = 0
    print('len(temp_seq):', len(temp_seq[k]),len(temp_seq[k+1]))
    embeddings = embedder.embed_many([str(s.seq) for s in temp_seq])
    embeddings = list(embeddings)
    print('embeddings[:1]:', embeddings[:1])
    for h in range(2):
        temp = pd.DataFrame(embeddings[h])
        print(j+h,sequences_names[j+h])
        file_name_path = '../data/MissenseM_ALL_ESM1bEmbedder/'+sequences_names[j+h]+'_ESM1bEmbedder.txt'
        temp.to_csv(file_name_path)
    print(file_name_path+' save sucess !')