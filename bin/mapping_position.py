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
            if len(i)>1:
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
    for t_pos in range(len(d_trim.keys())):
        for ori_pos in range(len(d_ori.keys())):
            if ori_pos>last_o and t_pos>last_t and d_ori[ori_pos]==d_trim[t_pos]:
                pos_map[t_pos]=ori_pos
                last_t=t_pos
                last_o=ori_pos
    return pos_map

def generating_3d_input_file(pos_map, fam,tr_alg):
    len_aln=len(pos_map.keys())
    out_file=str(fam)+'_'+str(tr_alg)+'_selected_columns.txt'
    out=open(out_file,'w')
    pairs=[]
    for p1 in range(len_aln):
        for p2 in range(p1+1, len_aln):
            #pairs.append((tmp_list[p1], tmp_list[p2]))
            pairs.append((int(pos_map[p1])+1,int(pos_map[p2])+1 ))
    for t in pairs:
        out.write(str(t[0])+' '+str(t[1])+'\n')
    out.close()


if __name__ == '__main__':
    msa_file=sys.argv[1]
    ori_msa=sys.argv[2]
    fam=sys.argv[3]
    tr_alg=sys.argv[4]
    lid1,trimmed_msa=import_msa(msa_file)
    lid2,original_msa=import_msa(ori_msa)
    pos_map=position_mapping(trimmed_msa, original_msa, lid1)
    generating_3d_input_file(pos_map, fam,tr_alg)

#    python mapping_position.py ${trimmed_msa} ${original_msa} ${id} trimal
