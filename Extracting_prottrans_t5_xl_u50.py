#!/usr/bin/env python
# coding: utf-8

def convert_three_letter_to_one_letter(three_letter_code):
    three_letter_code = three_letter_code.upper()
    amino_acids = {
        'ALA': 'A',
        'ARG': 'R',
        'ASN': 'N',
        'ASP': 'D',
        'CYS': 'C',
        'GLU': 'E',
        'GLN': 'Q',
        'GLY': 'G',
        'HIS': 'H',
        'ILE': 'I',
        'LEU': 'L',
        'LYS': 'K',
        'MET': 'M',
        'PHE': 'F',
        'PRO': 'P',
        'SER': 'S',
        'THR': 'T',
        'TRP': 'W',
        'TYR': 'Y',
        'VAL': 'V' }
    if three_letter_code in amino_acids:
        return amino_acids[three_letter_code]
    else:
        return 'Unknown'
print('define finish!')
a='MDKGRAAKVCHHADCQQLHRRGPLNLCEACDSKFHSTMHYDGHVRFDLPPQGSVLARNVSTRSCPPRTSPAVDLEEEEEESSVDGKGDRKSTGLKLSKKKARRRHTDDPSKECFTLKFDLNVDIETEIVPAMKKKSLGEVLLPVFERKGIALGKVDIYLDQSNTPLSLTFEAYRFGGHYLRVKAPAKPGDEGKVEQGMKDSKSLSLPILRPAGTGPPALERVDAQSRRESLDILAPGRRRKNMSEFLGEASIPGQEPPTPSSCSLPSGSSGSTNTGDSWKNRAASRFSGFFSSGPSTSAFGREVDKMEQLEGKLHTYSLFGLPRLPRGLRFDHDSWEEEYDEDEDEDNACLRLEDSWRELIDGHEKLTRRQCHQQEAVWELLHTEASYIRKLRVIINLFLCCLLNLQESGLLCEVEAERLFSNIPEIAQLHRRLWASVMAPVLEKARRTRALLQPGDFLKGFKMFGSLFKPYIRYCMEEEGCMEYMRGLLRDNDLFRAYITWAEKHPQCQRLKLSDMLAKPHQRLTKYPLLLKSVLRKTEEPRAKEAVVAMIGSVERFIHHVNACMRQRQERQRLAAVVSRIDAYEVVESSSDEVDKLLKEFLHLDLTAPIPGASPEETRQLLLEGSLRMKEGKDSKMDVYCFLFTDLLLVTKAVKKAERTRVIRPPLLVDKIVCRELRDPGSFLLIYLNEFHSAVGAYTFQASGQALCRGWVDTIYNAQNQLQQLRAQEPPGSQQPLQSLEEEEDEQEEEEEEEEEEEEGEDSGTSAASSPTIMRKSSGSPDSQHCASDGSTETLAMVVVEPGDTLSSPEFDSGPFSSQSDETSLSTTASSATPTSELLPLGPVDGRSCSMDSAYGTLSPTSLQDFVAPGPMAELVPRAPESPRVPSPPPSPRLRRRTPVQLLSCPPHLLKSKSEASLLQLLAGAGTHGTPSAPSRSLSELCLAVPAPGIRTQGSPQEAGPSWDCRGAPSPGSGPGLVGCLAGEPAGSHRKRCGDLPSGASPRVQPEPPPGVSAQHRKLTLAQLYRIRTTLLLNSTLTASEV'
a[731]
print(convert_three_letter_to_one_letter('ALA'))  # Output: A
print(convert_three_letter_to_one_letter('LEU'))  # Output: L
print(convert_three_letter_to_one_letter('XYZ'))  # Output: Unknown

def replace_amino_acid(sequence, wild_type, pos, mutant_AA):
    pos -= 1  # 转换为从零开始的索引
    sequence_list = list(sequence)
    if sequence_list[pos] == wild_type:
        sequence_list[pos] = mutant_AA
        mutated_sequence = ''.join(sequence_list)
        return mutated_sequence
    else:
        print("指定位置处的字母与野生型不匹配")
        return -1
print('define finish!')

import pandas as pd
import numpy as np
from pandas import DataFrame
from sklearn.preprocessing import StandardScaler, normalize
data = pd.read_excel('../data/5----all_MMdata_pro_seq_simple_label_var_2.xlsx', sheet_name='simple')
data.head(5)
Sequence = data['Sequence'].tolist()
ENSP = data['ENSP'].tolist()
Name = data['Name'].tolist()
Variant_AA_Change_original = data['Variant_AA_Change_original'].tolist()
Pos = data['Pos'].tolist()
Label = data['Label'].tolist()
len(Label)

mutantAA = []
wild_typeAA = []
for i in range(len(Variant_AA_Change_original )):
# for i in range(10):
    temp =Variant_AA_Change_original[i]
    tp2 = temp.split('.')[1]
    tp3 = tp2[:3]
    tp_wt_AA = convert_three_letter_to_one_letter(tp3)
    #print(tp2, tp3, tp_wt_AA)
    wild_typeAA.append(tp_wt_AA)
    tp = temp[-3:]
    tpAA = convert_three_letter_to_one_letter(tp)
    #print(tp, tpAA)
    mutantAA.append(tpAA)
print(mutantAA[:10], len(mutantAA),wild_typeAA[:10], len(wild_typeAA))
Seq = []
for i in range(len(Sequence)):
# for i in range(10):
    tp_seq = Sequence[i]
    tp_wtAA = wild_typeAA[i]
    tp_pos = Pos[i]-1
    tp_mutantAA = mutantAA[i]
    seq = replace_amino_acid(tp_seq, tp_wtAA, tp_pos, tp_mutantAA)
    Seq.append(seq)
print(len(Seq))

k = 0
for i in range(len(Seq)):
    temp = Seq[i]
    if(temp == -1):
        k = k+1
    else:
        continue
print(k)

from Bio import ExPASy
from Bio import SwissProt

protein_id = "ENSP00000368678"
with ExPASy.get_sprot_raw(protein_id) as handle:
    record = SwissProt.read(handle)
    sequence = record.sequence
print(sequence)
len_seq = []
for i in range(len(Sequence)):
    if(Sequence[i]=='-'):
        len_seq.append(0)
    else:
        len_seq.append(len(Sequence[i]))
len(len_seq)
SEQ = Sequence
print(len(ENSP),len(Sequence))

file = open('../data/00--independent test/3--VEP_output/5----Independent_MMdata_pro_seq_4_embeddings.txt','w')
for y in range(len(ENSP)):
    file.write('>'+ENSP[y]+'\n')
    file.write(SEQ[y]+'\n')
file.close()
print('save sucess !')

import os
path ="../data/MissenseM_Independent_ProtTransT5XLU50Embedder/"
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
for record in SeqIO.parse("../data/5----Independent_MMdata_pro_seq_4_embeddings.fasta", "fasta"):
    if (record.id not in processed_name):
        sequences.append(record)
        sequences_names.append(record.id)
    else:
        continue
    print('name--seq length:\t', record.id, len(record.seq))
print('length：', len(sequences),len(sequences_names))
embedder = ProtTransT5XLU50Embedder()
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
        file_name_path = '../data/MissenseM_Independent_ProtTransT5XLU50Embedder/'+sequences_names[j+h]+'_ProtTransT5XLU50Embedder.txt'
        temp.to_csv(file_name_path)
    print(file_name_path+' save sucess !')
import numpy as np
from bio_embeddings.embed import ProtTransT5XLU50Embedder
sequence = "MALLAMHSWRWAAAAAAFEKRRHSAILIRPLVSVSGSGPQWRPHQLGALGTARAYQQIPESLKSITWQRLGKGNSGQFLDAAKALQVWPLIEKRTCWHGHAGGGLHTDPKEGLKDVDTRKIIKAMLSYVWPKDRPDLRARVAISLGFLGGAKAMNIVVPFMFKYAVDSLNQMSGNMLNLSDAPNTVATMATAVLIGYGVSRAGAAFFNEVRNAVFGKVAQNSIRRIAKNVFLHLHNLDLGFHLSRQTGALSKAIDRGTRGISFVLSALVFNLLPIMFEVMLVSGVLYYKCGAQFALVTLGTLGTYTAFTVAVTRWRTRFRIEMNKADNDAGNAAIDSLLNYETVKYFNNERYEAQRYDGFLKTYETASLKSTSTLAMLNFGQSAIFSVGLTAIMVLASQGIVAGTLTVGDLVMVNGLLFQLSLPLNFLGTVYRETRQALIDMNTLFTLLKVDTQIKDKVMASPLQITPQTATVAFDNVHFEYIEGQKVLSGISFEVPAGKKVAIVGGSGSGKSTIVRLLFRFYEPQKGSIYLAGQNIQDVSLESLRRAVGVVPQDAVLFHNTIYYNLLYGNISASPEEVYAVAKLAGLHDAILRMPHGYDTQVGERGLKLSGGEKQRVAIARAILKDPPVILYDEATSSLDSITEETILGAMKDVVKHRTSIFIAHRLSTVVDADEIIVLDQGKVAERGTHHGLLANPHSIYSEMWHTQSSRVQNHDNPKWEAKKENISKEEERKKLQEEIVNSVKGCGNCSCMALLAMHSWRWAAAAAAFEKRRHSAILIRPLVSVSGSGPQWRPHQLGALGTARAYQQIPESLKSITWQRLGKGNSGQFLDAAKALQVWPLIEKRTCWHGHAGGGLHTDPKEGLKDVDTRKIIKAMLSYVWPKDRPDLRARVAISLGFLGGAKAMNIVVPFMFKYAVDSLNQMSGNMLNLSDAPNTVATMATAVLIGYGVSRAGAAFFNEVRNAVFGKVAQNSIRRIAKNVFLHLHNLDLGFHLSRQTGALSKAIDRGTRGISFVLSALVFNLLPIMFEVMLVSGVLYYKCGAQFALVTLGTLGTYTAFTVAVTRWRTRFRIEMNKADNDAGNAAIDSLLNYETVKYFNNERYEAQRYDGFLKTYETASLKSTSTLAMLNFGQSAIFSVGLTAIMVLASQGIVAGTLTVGDLVMVNGLLFQLSLPLNFLGTVYRETRQALIDMNTLFTLLKVDTQIKDKVMASPLQITPQTATVAFDNVHFEYIEGQKVLSGISFEVPAGKKVAIVGGSGSGKSTIVRLLFRFYEPQKGSIYLAGQNIQDVSLESLRRAVGVVPQDAVLFHNTIYYNLLYGNISASPEEVYAVAKLAGLHDAILRMPHGYDTQVGERGLKLSGGEKQRVAIARAILKDPPVILYDEATSSLDSITEETILGAMKDVVKHRTSIFIAHRLSTVVDADEIIVLDQGKVAERGTHHGLLANPHSIYSEMWHTQSSRVQNHDNPKWEAKKENISKEEERKKLQEEIVNSVKGCGNCSC"
print('sequence length:\t', len(sequence))
embedder = ProtTransT5XLU50Embedder()
print('embedder:\t', embedder)
embeddings = embedder.embed(sequence)
import pandas as pd
embeddings_save_1 = pd.DataFrame(embeddings)
embeddings_save_1.to_csv('../data/MissenseM_Independent_ProtTransT5XLU50Embedder/ENSP00000253577_ProtTransT5XLU50Embedder.txt')
print('file save sucess !')
import numpy as np
from bio_embeddings.embed import ProtTransT5XLU50Embedder
sequence = "MEELVVEVRGSNGAFYKAFVKDVHEDSITVAFENNWQPDRQIPFHDVRFPPPVGYNKDINESDEVEVYSRANEKEPCCWWLAKVRMIKGEFYVIEYAACDATYNEIVTIERLRSVNPNKPATKDTFHKIKLDVPEDLRQMCAKEAAHKDFKKAVGAFSVTYDPENYQLVILSINEVTSKRAHMLIDMHFRSLRTKLSLIMRNEEASKQLESSRQLASRFHEQFIVREDLMGLAIGTHGANIQQARKVPGVTAIDLDEDTCTFHIYGEDQDAVKKARSFLEFAEDVIQVPRNLVGKVIGKNGKLIQEIVDKSGVVRVRIEAENEKNVPQEEEIMPPNSLPSNNSRVGPNAPEEKKHLDIKENSTHFSQPNSTKVQRVLVASSVVAGESQKPELKAWQGMVPFVFVGTKDSIANATVLLDYHLNYLKEVDQLRLERLQIDEQLRQIGASSRPPPNRTDKEKSYVTDDGQGMGRGSRPYRNRGHGRRGPGYTSGTNSEASNASETESDHRDELSDWSLAPTEEERESFLRRGDGRRRGGGGRGQGGRGRGGGFKGNDDHSRTDNRPRNPREAKGRTTDGSLQIRVDCNNERSVHTKTLQNTSSEGSRLRTGKDRNQKKEKPDSVDGQQPLVNGVP"
print('sequence length:\t', len(sequence))
embedder = ProtTransT5XLU50Embedder()
print('embedder:\t', embedder)
embeddings = embedder.embed(sequence)
import pandas as pd
embeddings_save_1 = pd.DataFrame(embeddings)
embeddings_save_1.to_csv('../data/MissenseM_Independent_ProtTransT5XLU50Embedder/ENSP00000359506_ProtTransT5XLU50Embedder.txt')
print('file save sucess !')