#!/usr/bin/env python
# coding: utf-8

import pandas as pd
open_diff = open('./homo_sapiens_variation.txt/homo_sapiens_variation.txt','r',encoding='UTF-8')
diff_line = open_diff.readlines()
line_list = []
for line in diff_line:
    line_list.append(line)
count = len(line_list)

n = 5000000
print('hang number：', count)
diff_match_split = [line_list[i:i+n] for i in range(0, len(line_list),n)]

for i,j in zip(range(0, int(count/n+1)), range(0, int(count/n+1))):
    with open('./homo_sapiens_variation.txt/%d.txt'%j, 'w+',encoding='utf8') as temp:
        for line in diff_match_split[i]:
            temp.write(line)
print('the file number after spliting：', i+1)
open_diff.close()


missense = []
f0 = open('./homo_sapiens_variation.txt/0.txt','r',encoding='UTF-8')
f0_lines = f0.readlines()
for i in range(len(f0_lines)):
    temp = f0_lines[i]
    temp = temp.rstrip('\n')
    if("missense variant" in temp or "missense" in temp):
        if ("Likely benign"in temp or "Benign" in temp or "Pathogenic" in temp or "Likely pathogenic" in temp):
            missense.append(temp)
    else:
        continue
print(len(missense))
f0.close()

f1 = open('./homo_sapiens_variation.txt/1.txt','r',encoding='UTF-8')
f1_lines = f1.readlines()
for i in range(len(f1_lines)):
    temp = f1_lines[i]
    temp = temp.rstrip('\n')
    if("missense variant" in temp or "missense" in temp):
        if ("Likely benign"in temp or "Benign" in temp or "Pathogenic" in temp or "Likely pathogenic" in temp):
            missense.append(temp)
    else:
        continue
print(len(missense))
f1.close()

f2 = open('./homo_sapiens_variation.txt/2.txt','r',encoding='UTF-8')
f2_lines = f2.readlines()
for i in range(len(f2_lines)):
    temp = f2_lines[i]
    temp = temp.rstrip('\n')
    if("missense variant" in temp or "missense" in temp):
        if ("Likely benign"in temp or "Benign" in temp or "Pathogenic" in temp or "Likely pathogenic" in temp):
            missense.append(temp)
    else:
        continue
print(len(missense))
f2.close()

f3 = open('./homo_sapiens_variation.txt/3.txt','r',encoding='UTF-8')
f3_lines = f3.readlines()
for i in range(len(f3_lines)):
    temp = f3_lines[i]
    temp = temp.rstrip('\n')
    if("missense variant" in temp or "missense" in temp):
        if ("Likely benign"in temp or "Benign" in temp or "Pathogenic" in temp or "Likely pathogenic" in temp):
            missense.append(temp)
    else:
        continue
print(len(missense))
f3.close()

f4 = open('./homo_sapiens_variation.txt/4.txt','r',encoding='UTF-8')
f4_lines = f4.readlines()
for i in range(len(f4_lines)):
    temp = f4_lines[i]
    temp = temp.rstrip('\n')
    if("missense variant" in temp or "missense" in temp):
        if ("Likely benign"in temp or "Benign" in temp or "Pathogenic" in temp or "Likely pathogenic" in temp):
            missense.append(temp)
    else:
        continue
print(len(missense))
f4.close()

f5 = open('./homo_sapiens_variation.txt/5.txt','r',encoding='UTF-8')
f5_lines = f5.readlines()
for i in range(len(f5_lines)):
    temp = f5_lines[i]
    temp = temp.rstrip('\n')
    if("missense variant" in temp or "missense" in temp):
        if ("Likely benign"in temp or "Benign" in temp or "Pathogenic" in temp or "Likely pathogenic" in temp):
            missense.append(temp)
    else:
        continue
print(len(missense))
f5.close()

f6 = open('./homo_sapiens_variation.txt/6.txt','r',encoding='UTF-8')
f6_lines = f6.readlines()
for i in range(len(f6_lines)):
    temp = f6_lines[i]
    temp = temp.rstrip('\n')
    if("missense variant" in temp or "missense" in temp):
        if ("Likely benign"in temp or "Benign" in temp or "Pathogenic" in temp or "Likely pathogenic" in temp):
            missense.append(temp)
    else:
        continue
print(len(missense))
f6.close()

f7 = open('./homo_sapiens_variation.txt/7.txt','r',encoding='UTF-8')
f7_lines = f7.readlines()
for i in range(len(f7_lines)):
    temp = f7_lines[i]
    temp = temp.rstrip('\n')
    if("missense variant" in temp or "missense" in temp):
        if ("Likely benign"in temp or "Benign" in temp or "Pathogenic" in temp or "Likely pathogenic" in temp):
            missense.append(temp)
    else:
        continue
print(len(missense))
f7.close()
print('finish !')


import pandas as pd
df_missense = pd.DataFrame(missense)
df_missense.to_csv('./homo_sapiens_variation.txt/0__missense_benign_pathogenic_697797.txt',index='False')
print('save sucessfully !')



import pandas as pd
data = pd.read_excel('./0--homo_sapiens_variation.txt/0__missense_benign_pathogenic_697797.xlsx')
data.shape
Clinical_Significance_1 = data['Clinical Significance_1'].tolist()
Clinical_Significance_2 = data['Clinical Significance_2'].tolist()
Clinical_Significance_3 = data['Clinical Significance_3'].tolist()
Clinical_Significance_4 = data['Clinical Significance_4'].tolist()
Clinical_Significance_5 = data['Clinical Significance_5'].tolist()
Clinical_Significance_6 = data['Clinical Significance_6'].tolist()
print(len(Clinical_Significance_1), len(Clinical_Significance_6))
Clinical_Significance = []
for i in range(len(Clinical_Significance_1)):
    temp = str(Clinical_Significance_1[i])+str(Clinical_Significance_2[i])+str(Clinical_Significance_3[i])+str(Clinical_Significance_4[i])+str(Clinical_Significance_5[i])+str(Clinical_Significance_6[i])
    #print(temp)
    if ('Pathogenic' in temp or 'Likely pathogenic' in temp):
        Clinical_Significance.append('Pathogenic')
    elif ('Benign' in temp or 'Likely benign' in temp):
        Clinical_Significance.append('Benign')
    else:
        Clinical_Significance.append('-')
        print(temp) 
print('len(Clinical_Significance):\t', len(Clinical_Significance))

data_save = data
data_save['Clinical_Significance'] = pd.DataFrame(Clinical_Significance)
data_save.to_excel('./0--622270_original_MM.xlsx', index=False)
print('save sucess !')