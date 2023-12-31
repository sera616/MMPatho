#!/usr/bin/env python
# coding: utf-8

# In[6]:
import pandas as pd
data  = pd.read_excel('./0--622270_original_MM.xlsx')
# data.head(5)
print(data.shape)
print(data.columns)
Chromosome_Coordinate = data['Chromosome Coordinate'].tolist()
print('len(Chromosome_Coordinate):\t', len(Chromosome_Coordinate))

print('----- 1_5000 ----')
vcf_1_5000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/1-5000.vcf')
print(len(vcf_1_5000), type(vcf_1_5000))
print(vcf_1_5000.columns)
CHROM_1_5000 = vcf_1_5000['#CHROM'].tolist()
POS_1_5000 = vcf_1_5000['POS'].tolist()
ID_1_5000 = vcf_1_5000['ID'].tolist()
REF_1_5000 = vcf_1_5000['REF'].tolist()
ALT_1_5000 = vcf_1_5000['ALT'].tolist()
QUAL_1_5000 = vcf_1_5000['QUAL'].tolist()
FILTER_1_5000 = vcf_1_5000['FILTER'].tolist()
INFO_1_5000 = vcf_1_5000['INFO'].tolist()
print(len(CHROM_1_5000),len(POS_1_5000),len(ID_1_5000),len(REF_1_5000),len(ALT_1_5000))
print(len(QUAL_1_5000),len(FILTER_1_5000),len(INFO_1_5000))

print('\n----- 5001_10000 ----')
vcf_5001_10000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/5001-10000.vcf')
print(len(vcf_5001_10000), type(vcf_5001_10000))
print(vcf_5001_10000.columns)
CHROM_5001_10000 = vcf_5001_10000['#CHROM'].tolist()
POS_5001_10000 = vcf_5001_10000['POS'].tolist()
ID_5001_10000 = vcf_5001_10000['ID'].tolist()
REF_5001_10000 = vcf_5001_10000['REF'].tolist()
ALT_5001_10000 = vcf_5001_10000['ALT'].tolist()
QUAL_5001_10000 = vcf_5001_10000['QUAL'].tolist()
FILTER_5001_10000 = vcf_5001_10000['FILTER'].tolist()
INFO_5001_10000 = vcf_5001_10000['INFO'].tolist()
print(len(CHROM_5001_10000),len(POS_5001_10000),len(ID_5001_10000),len(REF_5001_10000),len(ALT_5001_10000))
print(len(QUAL_5001_10000),len(FILTER_5001_10000),len(INFO_5001_10000))

print('\n----- 10001_15000 ----')
vcf_10001_15000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/10001-15000.vcf')
print(len(vcf_10001_15000), type(vcf_10001_15000))
print(vcf_10001_15000.columns)
CHROM_10001_15000 = vcf_10001_15000['#CHROM'].tolist()
POS_10001_15000 = vcf_10001_15000['POS'].tolist()
ID_10001_15000 = vcf_10001_15000['ID'].tolist()
REF_10001_15000 = vcf_10001_15000['REF'].tolist()
ALT_10001_15000 = vcf_10001_15000['ALT'].tolist()
QUAL_10001_15000 = vcf_10001_15000['QUAL'].tolist()
FILTER_10001_15000 = vcf_10001_15000['FILTER'].tolist()
INFO_10001_15000 = vcf_10001_15000['INFO'].tolist()
print(len(CHROM_10001_15000),len(POS_10001_15000),len(ID_10001_15000),len(REF_10001_15000),len(ALT_10001_15000))
print(len(QUAL_10001_15000),len(FILTER_10001_15000),len(INFO_10001_15000))

print('\n----- 15001_20000 ----')
vcf_15001_20000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/15001-20000.vcf')
print(len(vcf_15001_20000), type(vcf_15001_20000))
print(vcf_15001_20000.columns)
CHROM_15001_20000 = vcf_15001_20000['#CHROM'].tolist()
POS_15001_20000 = vcf_15001_20000['POS'].tolist()
ID_15001_20000 = vcf_15001_20000['ID'].tolist()
REF_15001_20000 = vcf_15001_20000['REF'].tolist()
ALT_15001_20000 = vcf_15001_20000['ALT'].tolist()
QUAL_15001_20000 = vcf_15001_20000['QUAL'].tolist()
FILTER_15001_20000 = vcf_15001_20000['FILTER'].tolist()
INFO_15001_20000 = vcf_15001_20000['INFO'].tolist()
print(len(CHROM_15001_20000),len(POS_15001_20000),len(ID_15001_20000),len(REF_15001_20000),len(ALT_15001_20000))
print(len(QUAL_15001_20000),len(FILTER_15001_20000),len(INFO_15001_20000))

print('\n----- 20001_25000 ----')
vcf_20001_25000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/20001-25000.vcf')
print(len(vcf_20001_25000), type(vcf_20001_25000))
print(vcf_20001_25000.columns)
CHROM_20001_25000 = vcf_20001_25000['#CHROM'].tolist()
POS_20001_25000 = vcf_20001_25000['POS'].tolist()
ID_20001_25000 = vcf_20001_25000['ID'].tolist()
REF_20001_25000 = vcf_20001_25000['REF'].tolist()
ALT_20001_25000 = vcf_20001_25000['ALT'].tolist()
QUAL_20001_25000 = vcf_20001_25000['QUAL'].tolist()
FILTER_20001_25000 = vcf_20001_25000['FILTER'].tolist()
INFO_20001_25000 = vcf_20001_25000['INFO'].tolist()
print(len(CHROM_20001_25000),len(POS_20001_25000),len(ID_20001_25000),len(REF_20001_25000),len(ALT_20001_25000))
print(len(QUAL_20001_25000),len(FILTER_20001_25000),len(INFO_20001_25000))

print('\n----- 25001_30000 ----')
vcf_25001_30000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/25001-30000.vcf')
print(len(vcf_25001_30000), type(vcf_25001_30000))
print(vcf_25001_30000.columns)
CHROM_25001_30000 = vcf_25001_30000['#CHROM'].tolist()
POS_25001_30000 = vcf_25001_30000['POS'].tolist()
ID_25001_30000 = vcf_25001_30000['ID'].tolist()
REF_25001_30000 = vcf_25001_30000['REF'].tolist()
ALT_25001_30000 = vcf_25001_30000['ALT'].tolist()
QUAL_25001_30000 = vcf_25001_30000['QUAL'].tolist()
FILTER_25001_30000 = vcf_25001_30000['FILTER'].tolist()
INFO_25001_30000 = vcf_25001_30000['INFO'].tolist()
print(len(CHROM_25001_30000),len(POS_25001_30000),len(ID_25001_30000),len(REF_25001_30000),len(ALT_25001_30000))
print(len(QUAL_25001_30000),len(FILTER_25001_30000),len(INFO_25001_30000))


print('\n----- 30001_35000 ----')
vcf_30001_35000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/30001-35000.vcf')
print(len(vcf_30001_35000), type(vcf_30001_35000))
print(vcf_30001_35000.columns)
CHROM_30001_35000 = vcf_30001_35000['#CHROM'].tolist()
POS_30001_35000 = vcf_30001_35000['POS'].tolist()
ID_30001_35000 = vcf_30001_35000['ID'].tolist()
REF_30001_35000 = vcf_30001_35000['REF'].tolist()
ALT_30001_35000 = vcf_30001_35000['ALT'].tolist()
QUAL_30001_35000 = vcf_30001_35000['QUAL'].tolist()
FILTER_30001_35000 = vcf_30001_35000['FILTER'].tolist()
INFO_30001_35000 = vcf_30001_35000['INFO'].tolist()
print(len(CHROM_30001_35000),len(POS_30001_35000),len(ID_30001_35000),len(REF_30001_35000),len(ALT_30001_35000))
print(len(QUAL_30001_35000),len(FILTER_30001_35000),len(INFO_30001_35000))

print('\n----- 35001_40000 ----')
vcf_35001_40000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/35001-40000.vcf')
print(len(vcf_35001_40000), type(vcf_35001_40000))
print(vcf_35001_40000.columns)
CHROM_35001_40000 = vcf_35001_40000['#CHROM'].tolist()
POS_35001_40000 = vcf_35001_40000['POS'].tolist()
ID_35001_40000 = vcf_35001_40000['ID'].tolist()
REF_35001_40000 = vcf_35001_40000['REF'].tolist()
ALT_35001_40000 = vcf_35001_40000['ALT'].tolist()
QUAL_35001_40000 = vcf_35001_40000['QUAL'].tolist()
FILTER_35001_40000 = vcf_35001_40000['FILTER'].tolist()
INFO_35001_40000 = vcf_35001_40000['INFO'].tolist()
print(len(CHROM_35001_40000),len(POS_35001_40000),len(ID_35001_40000),len(REF_35001_40000),len(ALT_35001_40000))
print(len(QUAL_35001_40000),len(FILTER_35001_40000),len(INFO_35001_40000))

print('\n----- 40001_45000 ----')
vcf_40001_45000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/40001-45000.vcf')
print(len(vcf_40001_45000), type(vcf_40001_45000))
print(vcf_40001_45000.columns)
CHROM_40001_45000 = vcf_40001_45000['#CHROM'].tolist()
POS_40001_45000 = vcf_40001_45000['POS'].tolist()
ID_40001_45000 = vcf_40001_45000['ID'].tolist()
REF_40001_45000 = vcf_40001_45000['REF'].tolist()
ALT_40001_45000 = vcf_40001_45000['ALT'].tolist()
QUAL_40001_45000 = vcf_40001_45000['QUAL'].tolist()
FILTER_40001_45000 = vcf_40001_45000['FILTER'].tolist()
INFO_40001_45000 = vcf_40001_45000['INFO'].tolist()
print(len(CHROM_40001_45000),len(POS_40001_45000),len(ID_40001_45000),len(REF_40001_45000),len(ALT_40001_45000))
print(len(QUAL_40001_45000),len(FILTER_40001_45000),len(INFO_40001_45000))

print('\n----- 45001_50000 ----')
vcf_45001_50000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/45001-50000.vcf')
print(len(vcf_45001_50000), type(vcf_45001_50000))
print(vcf_45001_50000.columns)
CHROM_45001_50000 = vcf_45001_50000['#CHROM'].tolist()
POS_45001_50000 = vcf_45001_50000['POS'].tolist()
ID_45001_50000 = vcf_45001_50000['ID'].tolist()
REF_45001_50000 = vcf_45001_50000['REF'].tolist()
ALT_45001_50000 = vcf_45001_50000['ALT'].tolist()
QUAL_45001_50000 = vcf_45001_50000['QUAL'].tolist()
FILTER_45001_50000 = vcf_45001_50000['FILTER'].tolist()
INFO_45001_50000 = vcf_45001_50000['INFO'].tolist()
print(len(CHROM_45001_50000),len(POS_45001_50000),len(ID_45001_50000),len(REF_45001_50000),len(ALT_45001_50000))
print(len(QUAL_45001_50000),len(FILTER_45001_50000),len(INFO_45001_50000))

print('\n----- 50001_55000 ----')
vcf_50001_55000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/50001-55000.vcf')
print(len(vcf_50001_55000), type(vcf_50001_55000))
print(vcf_50001_55000.columns)
CHROM_50001_55000 = vcf_50001_55000['#CHROM'].tolist()
POS_50001_55000 = vcf_50001_55000['POS'].tolist()
ID_50001_55000 = vcf_50001_55000['ID'].tolist()
REF_50001_55000 = vcf_50001_55000['REF'].tolist()
ALT_50001_55000 = vcf_50001_55000['ALT'].tolist()
QUAL_50001_55000 = vcf_50001_55000['QUAL'].tolist()
FILTER_50001_55000 = vcf_50001_55000['FILTER'].tolist()
INFO_50001_55000 = vcf_50001_55000['INFO'].tolist()
print(len(CHROM_50001_55000),len(POS_50001_55000),len(ID_50001_55000),len(REF_50001_55000),len(ALT_50001_55000))
print(len(QUAL_50001_55000),len(FILTER_50001_55000),len(INFO_50001_55000))

print('\n----- 55001_60000 ----')
vcf_55001_60000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/55001-60000.vcf')
print(len(vcf_55001_60000), type(vcf_55001_60000))
print(vcf_55001_60000.columns)
CHROM_55001_60000 = vcf_55001_60000['#CHROM'].tolist()
POS_55001_60000 = vcf_55001_60000['POS'].tolist()
ID_55001_60000 = vcf_55001_60000['ID'].tolist()
REF_55001_60000 = vcf_55001_60000['REF'].tolist()
ALT_55001_60000 = vcf_55001_60000['ALT'].tolist()
QUAL_55001_60000 = vcf_55001_60000['QUAL'].tolist()
FILTER_55001_60000 = vcf_55001_60000['FILTER'].tolist()
INFO_55001_60000 = vcf_55001_60000['INFO'].tolist()
print(len(CHROM_55001_60000),len(POS_55001_60000),len(ID_55001_60000),len(REF_55001_60000),len(ALT_55001_60000))
print(len(QUAL_55001_60000),len(FILTER_55001_60000),len(INFO_55001_60000))

print('\n----- 60001_65000 ----')
vcf_60001_65000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/60001-65000.vcf')
print(len(vcf_60001_65000), type(vcf_60001_65000))
print(vcf_60001_65000.columns)
CHROM_60001_65000 = vcf_60001_65000['#CHROM'].tolist()
POS_60001_65000 = vcf_60001_65000['POS'].tolist()
ID_60001_65000 = vcf_60001_65000['ID'].tolist()
REF_60001_65000 = vcf_60001_65000['REF'].tolist()
ALT_60001_65000 = vcf_60001_65000['ALT'].tolist()
QUAL_60001_65000 = vcf_60001_65000['QUAL'].tolist()
FILTER_60001_65000 = vcf_60001_65000['FILTER'].tolist()
INFO_60001_65000 = vcf_60001_65000['INFO'].tolist()
print(len(CHROM_60001_65000),len(POS_60001_65000),len(ID_60001_65000),len(REF_60001_65000),len(ALT_60001_65000))
print(len(QUAL_60001_65000),len(FILTER_60001_65000),len(INFO_60001_65000))

print('\n----- 65001_70000 ----')
vcf_65001_70000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/65001-70000.vcf')
print(len(vcf_65001_70000), type(vcf_65001_70000))
print(vcf_65001_70000.columns)
CHROM_65001_70000 = vcf_65001_70000['#CHROM'].tolist()
POS_65001_70000 = vcf_65001_70000['POS'].tolist()
ID_65001_70000 = vcf_65001_70000['ID'].tolist()
REF_65001_70000 = vcf_65001_70000['REF'].tolist()
ALT_65001_70000 = vcf_65001_70000['ALT'].tolist()
QUAL_65001_70000 = vcf_65001_70000['QUAL'].tolist()
FILTER_65001_70000 = vcf_65001_70000['FILTER'].tolist()
INFO_65001_70000 = vcf_65001_70000['INFO'].tolist()
print(len(CHROM_65001_70000),len(POS_65001_70000),len(ID_65001_70000),len(REF_65001_70000),len(ALT_65001_70000))
print(len(QUAL_65001_70000),len(FILTER_65001_70000),len(INFO_65001_70000))

print('\n----- 70001_75000 ----')
vcf_70001_75000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/70001-75000.vcf')
print(len(vcf_70001_75000), type(vcf_70001_75000))
print(vcf_70001_75000.columns)
CHROM_70001_75000 = vcf_70001_75000['#CHROM'].tolist()
POS_70001_75000 = vcf_70001_75000['POS'].tolist()
ID_70001_75000 = vcf_70001_75000['ID'].tolist()
REF_70001_75000 = vcf_70001_75000['REF'].tolist()
ALT_70001_75000 = vcf_70001_75000['ALT'].tolist()
QUAL_70001_75000 = vcf_70001_75000['QUAL'].tolist()
FILTER_70001_75000 = vcf_70001_75000['FILTER'].tolist()
INFO_70001_75000 = vcf_70001_75000['INFO'].tolist()
print(len(CHROM_70001_75000),len(POS_70001_75000),len(ID_70001_75000),len(REF_70001_75000),len(ALT_70001_75000))
print(len(QUAL_70001_75000),len(FILTER_70001_75000),len(INFO_70001_75000))

print('\n----- 75001_80000 ----')
vcf_75001_80000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/75001-80000.vcf')
print(len(vcf_75001_80000), type(vcf_75001_80000))
print(vcf_75001_80000.columns)
CHROM_75001_80000 = vcf_75001_80000['#CHROM'].tolist()
POS_75001_80000 = vcf_75001_80000['POS'].tolist()
ID_75001_80000 = vcf_75001_80000['ID'].tolist()
REF_75001_80000 = vcf_75001_80000['REF'].tolist()
ALT_75001_80000 = vcf_75001_80000['ALT'].tolist()
QUAL_75001_80000 = vcf_75001_80000['QUAL'].tolist()
FILTER_75001_80000 = vcf_75001_80000['FILTER'].tolist()
INFO_75001_80000 = vcf_75001_80000['INFO'].tolist()
print(len(CHROM_75001_80000),len(POS_75001_80000),len(ID_75001_80000),len(REF_75001_80000),len(ALT_75001_80000))
print(len(QUAL_75001_80000),len(FILTER_75001_80000),len(INFO_75001_80000))

print('\n----- 80001_85000 ----')
vcf_80001_85000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/80001-85000.vcf')
print(len(vcf_80001_85000), type(vcf_80001_85000))
print(vcf_80001_85000.columns)
CHROM_80001_85000 = vcf_80001_85000['#CHROM'].tolist()
POS_80001_85000 = vcf_80001_85000['POS'].tolist()
ID_80001_85000 = vcf_80001_85000['ID'].tolist()
REF_80001_85000 = vcf_80001_85000['REF'].tolist()
ALT_80001_85000 = vcf_80001_85000['ALT'].tolist()
QUAL_80001_85000 = vcf_80001_85000['QUAL'].tolist()
FILTER_80001_85000 = vcf_80001_85000['FILTER'].tolist()
INFO_80001_85000 = vcf_80001_85000['INFO'].tolist()
print(len(CHROM_80001_85000),len(POS_80001_85000),len(ID_80001_85000),len(REF_80001_85000),len(ALT_80001_85000))
print(len(QUAL_80001_85000),len(FILTER_80001_85000),len(INFO_80001_85000))

print('\n----- 85001_90000 ----')
vcf_85001_90000 = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/85001-90000.vcf')
print(len(vcf_85001_90000), type(vcf_85001_90000))
print(vcf_85001_90000.columns)
CHROM_85001_90000 = vcf_85001_90000['#CHROM'].tolist()
POS_85001_90000 = vcf_85001_90000['POS'].tolist()
ID_85001_90000 = vcf_85001_90000['ID'].tolist()
REF_85001_90000 = vcf_85001_90000['REF'].tolist()
ALT_85001_90000 = vcf_85001_90000['ALT'].tolist()
QUAL_85001_90000 = vcf_85001_90000['QUAL'].tolist()
FILTER_85001_90000 = vcf_85001_90000['FILTER'].tolist()
INFO_85001_90000 = vcf_85001_90000['INFO'].tolist()
print(len(CHROM_85001_90000),len(POS_85001_90000),len(ID_85001_90000),len(REF_85001_90000),len(ALT_85001_90000))
print(len(QUAL_85001_90000),len(FILTER_85001_90000),len(INFO_85001_90000))

print('\n----- 90001_final ----')
vcf_90001_final = pd.read_table('./1--Ensembl Variant Recorder__input_output/output/90001-final.vcf')
print(len(vcf_90001_final), type(vcf_90001_final))
print(vcf_90001_final.columns)
CHROM_90001_final = vcf_90001_final['#CHROM'].tolist()
POS_90001_final = vcf_90001_final['POS'].tolist()
ID_90001_final = vcf_90001_final['ID'].tolist()
REF_90001_final = vcf_90001_final['REF'].tolist()
ALT_90001_final = vcf_90001_final['ALT'].tolist()
QUAL_90001_final = vcf_90001_final['QUAL'].tolist()
FILTER_90001_final = vcf_90001_final['FILTER'].tolist()
INFO_90001_final = vcf_90001_final['INFO'].tolist()
print(len(CHROM_90001_final),len(POS_90001_final),len(ID_90001_final),len(REF_90001_final),len(ALT_90001_final))
print(len(QUAL_90001_final),len(FILTER_90001_final),len(INFO_90001_final))

CHROM = CHROM_1_5000 + CHROM_5001_10000 + CHROM_10001_15000 + CHROM_15001_20000 + CHROM_20001_25000 + CHROM_25001_30000 + CHROM_30001_35000 + CHROM_35001_40000 + CHROM_40001_45000 + CHROM_45001_50000 + CHROM_50001_55000 + CHROM_55001_60000 + CHROM_60001_65000 + CHROM_65001_70000 + CHROM_70001_75000 + CHROM_75001_80000 + CHROM_80001_85000 + CHROM_85001_90000 + CHROM_90001_final
POS = POS_1_5000 + POS_5001_10000 + POS_10001_15000 + POS_15001_20000 + POS_20001_25000 + POS_25001_30000 + POS_30001_35000 + POS_35001_40000 + POS_40001_45000 + POS_45001_50000 + POS_50001_55000 + POS_55001_60000 + POS_60001_65000 + POS_65001_70000 + POS_70001_75000 + POS_75001_80000 + POS_80001_85000 + POS_85001_90000 + POS_90001_final
ID = ID_1_5000 + ID_5001_10000 + ID_10001_15000 + ID_15001_20000 + ID_20001_25000 + ID_25001_30000 + ID_30001_35000 + ID_35001_40000 + ID_40001_45000 + ID_45001_50000 + ID_50001_55000 + ID_55001_60000 + ID_60001_65000 + ID_65001_70000 + ID_70001_75000 + ID_75001_80000 + ID_80001_85000 + ID_85001_90000 + ID_90001_final
REF = REF_1_5000 + REF_5001_10000 + REF_10001_15000 + REF_15001_20000 + REF_20001_25000 + REF_25001_30000 + REF_30001_35000 + REF_35001_40000 + REF_40001_45000 + REF_45001_50000 + REF_50001_55000 + REF_55001_60000 + REF_60001_65000 + REF_65001_70000 + REF_70001_75000 + REF_75001_80000 + REF_80001_85000 + REF_85001_90000 + REF_90001_final
ALT = ALT_1_5000 + ALT_5001_10000 + ALT_10001_15000 + ALT_15001_20000 + ALT_20001_25000 + ALT_25001_30000 + ALT_30001_35000 + ALT_35001_40000 + REF_40001_45000 + ALT_45001_50000 + ALT_50001_55000 + ALT_55001_60000 + ALT_60001_65000 + ALT_65001_70000 + ALT_70001_75000 + ALT_75001_80000 + ALT_80001_85000 + ALT_85001_90000 + ALT_90001_final
QUAL = QUAL_1_5000 + QUAL_5001_10000 + QUAL_10001_15000 + QUAL_15001_20000 + QUAL_20001_25000 + QUAL_25001_30000 + QUAL_30001_35000 + QUAL_35001_40000 + QUAL_40001_45000 + QUAL_45001_50000 + QUAL_50001_55000 + QUAL_55001_60000 + QUAL_60001_65000 + QUAL_65001_70000 + QUAL_70001_75000 + QUAL_75001_80000 + QUAL_80001_85000 + QUAL_85001_90000 + QUAL_90001_final
FILTER = FILTER_1_5000 +  FILTER_5001_10000 + FILTER_10001_15000 + FILTER_15001_20000 + FILTER_20001_25000 + FILTER_25001_30000 + FILTER_30001_35000 + FILTER_35001_40000 + FILTER_40001_45000 + FILTER_45001_50000 + FILTER_50001_55000 + FILTER_55001_60000 + FILTER_60001_65000 + FILTER_65001_70000 + FILTER_70001_75000 + FILTER_75001_80000 + FILTER_80001_85000 + FILTER_85001_90000 + FILTER_90001_final
INFO = INFO_1_5000 + INFO_5001_10000 + INFO_10001_15000 + INFO_15001_20000 + INFO_20001_25000 + INFO_25001_30000 + INFO_30001_35000 + INFO_35001_40000 + INFO_40001_45000 + INFO_45001_50000 + INFO_50001_55000 + INFO_55001_60000 + INFO_60001_65000 + INFO_65001_70000 + INFO_70001_75000 + INFO_75001_80000 + INFO_80001_85000 + INFO_85001_90000 + INFO_90001_final
print(len(CHROM),len(POS),len(ID),len(REF),len(ALT),len(FILTER),len(INFO))


def get_index_list(lst, target):
    a = [i for i,v in enumerate(lst) if v == target]
    return(a)

CHROM_keep = []
POS_keep = []
ID_keep = []
REF_keep = []
ALT_keep = []
QUAL_keep = []
FILTER_keep = []
INFO_keep = []
for i in range(len(Chromosome_Coordinate)):
#for i in range(100):
    temp = Chromosome_Coordinate[i]
    xiabiao_temp = []
    xiabiao_temp = get_index_list(ID, temp)
    print(xiabiao_temp,type(xiabiao_temp),len(xiabiao_temp))
    if(len(xiabiao_temp)>1):
        print('----have two index：\t')
    try:        
        xiabiao = ID.index(temp)
        CHROM_keep.append(CHROM[xiabiao])
        POS_keep.append(POS[xiabiao])
        ID_keep.append(ID[xiabiao])
        REF_keep.append(REF[xiabiao])
        ALT_keep.append(ALT[xiabiao])
        QUAL_keep.append(QUAL[xiabiao])
        FILTER_keep.append(FILTER[xiabiao])
        INFO_keep.append(INFO[xiabiao])     
    except ValueError:
        CHROM_keep.append('-')
        POS_keep.append('-')
        ID_keep.append('-')
        REF_keep.append('-')
        ALT_keep.append('-')
        QUAL_keep.append('-')
        FILTER_keep.append('-')
        INFO_keep.append('-')   

print('finish\n')
print(len(CHROM_keep),len(POS_keep),len(ID_keep),len(REF_keep),len(ALT_keep),len(FILTER_keep),len(INFO_keep))

data_save = data
data_save['CHROM_grch38'] = pd.DataFrame(CHROM_keep)
data_save['POS_grch38'] = pd.DataFrame(POS_keep)
data_save['ID_grch38'] = pd.DataFrame(ID_keep)
data_save['REF_grch38'] = pd.DataFrame(REF_keep)
data_save['ALT_grch38'] = pd.DataFrame(ALT_keep)
data_save['QUAL_grch38'] = pd.DataFrame(QUAL_keep)
data_save['FILTER_grch38'] = pd.DataFrame(FILTER_keep)
data_save['INFO_grch38'] = pd.DataFrame(INFO_keep)
data_save.to_excel('./1--91072_VEP_mapped_GRch38_MM.xlsx', index = False)
print('finish !')