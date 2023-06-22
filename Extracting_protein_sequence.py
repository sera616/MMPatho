#!/usr/bin/env python
# coding: utf-8

# In[1]:
import numpy as np
import pandas as pd
from pandas import DataFrame
from sklearn.preprocessing import StandardScaler, normalize
data = pd.read_csv('./4--Missense_Mutation_GRch37_for_VEP_info_output_mapped_features_temp4proseq.csv')
print(len(data))# 
data.head()
Uploaded_variation = data['Uploaded_variation'].tolist()
Location = data['Location'].tolist()
Allele = data['Allele'].tolist()
Gene = data['Gene'].tolist()
CHROM_grch37_original = data['CHROM_grch37_original'].tolist()
POS_grch37_original = data['POS_grch37_original'].tolist()
ID_grch37_original = data['ID_grch37_original'].tolist()
REF_grch37_original = data['REF_grch37_original'].tolist()
ALT_grch37_original = data['ALT_grch37_original'].tolist()
QUAL_grch37_original = data['QUAL_grch37_original'].tolist()
FILTER_grch37_original = data['FILTER_grch37_original'].tolist()
INFO_grch37_original = data['INFO_grch37_original'].tolist()
ENSP = data['ENSP'].tolist()
label = data['label'].tolist()
print(len(ENSP))

protein_seq = open('5----gall_MMdata_pro_seq.txt','w')
protein_seq.write('Uploaded_variation\t'+'Location\t'+'Allele\t'+'Gene\t'+
                  'CHROM_grch37_original\t'+'POS_grch37_original\t'+'ID_grch37_original\t'+'REF_grch37_original\t'+
                  'ALT_grch37_original\t'+'QUAL_grch37_original\t'+'FILTER_grch37_original\t'+'INFO_grch37_original\t'+
                  'ENSP\t'+'label\n')
import requests, sys
server = "https://rest.ensembl.org"
for i in range(len(ENSP)):    
    temp=ENSP[i]
    if (temp=='-'):
        protein_seq.write(str(Uploaded_variation[i])+'\t'+str(Location[i])+'\t'+str(Allele[i])+'\t'+str(Gene[i])+'\t'+
                  str(CHROM_grch37_original[i])+'\t'+str(POS_grch37_original[i])+'\t'+
                  str(ID_grch37_original[i])+'\t'+str(REF_grch37_original[i])+'\t'+str(ALT_grch37_original[i])+'\t'+
                  str(QUAL_grch37_original[i])+'\t'+str(FILTER_grch37_original[i])+'\t'+str(INFO_grch37_original[i])+'\t'+
                  str(ENSP[i])+'\t'+str(label[i])+'\t'+'-'+'\n')
        pass
    else:
        try:
            ext="/sequence/id/"+str(temp)+"?"
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
            if not r.ok:
                r.raise_for_status()
                sys.exit()
            decoded = r.json()
            print('-------'+temp+'---------')
            print('--repr(decoded)--',repr(decoded)) #type(repr(decoded)): str
            print('--decoded.sequence--',decoded['seq']) #type(decoded):dict    
            seq_temp = decoded['seq']
            protein_seq.write(str(Uploaded_variation[i])+'\t'+str(Location[i])+'\t'+str(Allele[i])+'\t'+str(Gene[i])+'\t'+
                  str(CHROM_grch37_original[i])+'\t'+str(POS_grch37_original[i])+'\t'+
                  str(ID_grch37_original[i])+'\t'+str(REF_grch37_original[i])+'\t'+str(ALT_grch37_original[i])+'\t'+
                  str(QUAL_grch37_original[i])+'\t'+str(FILTER_grch37_original[i])+'\t'+str(INFO_grch37_original[i])+'\t'+
                  str(ENSP[i])+'\t'+str(label[i])+'\t'+str(seq_temp)+'\n')
        except:
            protein_seq.write(str(Uploaded_variation[i])+'\t'+str(Location[i])+'\t'+str(Allele[i])+'\t'+str(Gene[i])+'\t'+
                              str(CHROM_grch37_original[i])+'\t'+str(POS_grch37_original[i])+'\t'+
                              str(ID_grch37_original[i])+'\t'+str(REF_grch37_original[i])+'\t'+str(ALT_grch37_original[i])+'\t'+
                              str(QUAL_grch37_original[i])+'\t'+str(FILTER_grch37_original[i])+'\t'+str(INFO_grch37_original[i])+'\t'+
                              str(ENSP[i])+'\t'+str(label[i])+'\t'+'-'+'\n')
            continue
protein_seq.close()
print('Extracting protein sequence finish!')