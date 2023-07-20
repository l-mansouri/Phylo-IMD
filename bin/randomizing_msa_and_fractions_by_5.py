#! usr/bin/env python
import sys
import random
import os

def import_msa(msa_file):
    '''importing the alignment'''
    msa={}
    lid=[]
    with open(msa_file) as f:
        for i in f:
            i=i.rstrip()
            if i[0]=='>':
                ID=i[1:]
                lid.append(ID)
                continue
            msa[ID]=msa.get(ID,'')+i
    return lid, msa

def position_mapping(trimmed_msa, original_msa, lid1):
    d_trim={}
    d_ori={}
    for i in range(len(trimmed_msa[lid1[0]])):
        s=''
        for j in range(len(lid1)):
            s=s+trimmed_msa[lid1[j]][i]
        d_trim[i]=s
    for i in range(len(original_msa[lid1[0]])):
        s=''
        for j in range(len(lid1)):
            s=s+original_msa[lid1[j]][i]
        d_ori[i]=s
    pos_map={}
    last_t=-1
    last_o=-1
    for t_pos in d_trim.keys():
        for ori_pos in d_ori.keys():
            if ori_pos>last_o and t_pos>last_t and d_ori[ori_pos]==d_trim[t_pos]:
                pos_map[t_pos]=ori_pos
                last_t=t_pos
                last_o=ori_pos
    return pos_map


def generating_random_msa(msa):
    '''randomizing the alignment'''
    #list=[0,..,len(aln)-1]
    #random.shuffle(list)
    randomized_list=[]
    random_msa=[]
    len_aln=len(msa[msa.keys()[0]])
    seeds=[50,682176,3473,2981,12,2445,9873,98457,497,9563]
    for rep in range(0,10):
        random.seed(seeds[rep])
        lis=[i for i in range(len_aln)]
        random.shuffle(lis)#list of shuffled positions
        randomized_list.append(lis[0:200])
        d={}
        for k in msa.iterkeys():
            #print k
            shuffled_seq=''
            for ind in lis[0:200]:
                #print ind
                shuffled_seq=shuffled_seq+msa[k][ind]#shuffled sequence for id k
                #print shuffled_seq
            d[k]=shuffled_seq
        random_msa.append(d)#list of dictionaries, each containing a shuffled alignment
        #print d
    return randomized_list, random_msa


def generating_1d_fractions(random_msa, fam, aln_type):
    #select fractions i+10 and write each one in a file
    len_aln=len(random_msa[0][random_msa[0].keys()[0]])
    for n in range(len(random_msa)):
        for i in range(5,len_aln,5):
            init=0
            fin=i+5
            out_file=str(fam)+'_'+str(aln_type)+'_random_msa_replicate.'+str(n)+'_with_'+str(fin)+'_columns.fa'
            out=open(out_file,'w')
            for k in random_msa[n].iterkeys():
                out.write('>'+str(k)+'\n')
                out.write(str(random_msa[n][k][init:fin])+'\n')
            out.close()

def generating_3d_input_file(randomized_list,pos_map, fam, aln_type):
    len_aln=len(randomized_list[0])
    results=[]
    for n in range(len(randomized_list)):
        lis=randomized_list[n]
        sel_pairs=[]
        for i in range(5, len_aln,5):
            fin=i+5
            tmp_list=lis[0:fin]
            out_file=str(fam)+'_'+str(aln_type)+'.random_column_pairs_replicate.'+str(n)+'_with_'+str(fin)+'_columns.txt'
            out=open(out_file,'w')
            pairs=[]
            tr_pairs=[]
            for p1 in range(len(tmp_list)):
                for p2 in range(p1+1, len(tmp_list)):
                    #pairs.append((tmp_list[p1], tmp_list[p2]))
                    tr_pairs.append((int(pos_map[tmp_list[p1]])+1,int(pos_map[tmp_list[p2]])+1 ))
            for t in tr_pairs:
                if t[0]<t[1]:
                    out.write(str(t[1])+' '+str(t[0])+'\n')
                else:
                    out.write(str(t[0])+' '+str(t[1])+'\n')
            out.close()



if __name__ == '__main__':
    msa_file=sys.argv[1]
    ori_msa=sys.argv[2]
    fam=sys.argv[3]
    aln_type=sys.argv[4]
    lid1,trimmed_msa=import_msa(msa_file)
    lid2,original_msa=import_msa(ori_msa)
    pos_map=position_mapping(trimmed_msa, original_msa, lid1)
    randomized_list, random_msa=generating_random_msa(trimmed_msa)
    generating_1d_fractions(random_msa, fam, aln_type)
    generating_3d_input_file(randomized_list, pos_map, fam, aln_type)

