# !/usr/bin/env python
#INPUT:
#SPLIT1,prova3dtrees.nwk,2,100,100000001000,original,1d
#
#OUTPUT:
#br_id ML_bs ME_bs IMD_bs in_ref_br
import sys

def import_split(splitfile):
    D={}
    for line in open(splitfile):
        line=line.rstrip().split(',')
        D[line[4]]=line[3]
    return D
#D={br_id: (bs)}

def comparing_references(rIMD, rME, rML):
    ref_branches=[]
    for k in rIMD.keys():
        if k in rME.keys() and k in rML.keys():
            ref_branches.append(k)
    return ref_branches


def comparing_3d_ME_ML(IMD, ME, ML, ref_branches, outfile):
    out=open(outfile, 'w')
    for k in IMD.keys():
        out.write(str(k)+' '+str(IMD[k])+' ')
        if k in ML.keys():
            out.write(str(ML[k])+' ')
        elif k not in ML.keys():
            out.write('0 ')
        if k in ME.keys():
            out.write(str(ME[k])+' ')
        elif k not in ME.keys():
            out.write('0 ')
        if k in ref_branches:
            out.write('1 \n')
        elif k not in ref_branches:
            out.write('0 \n')
    out.close()

if __name__ == '__main__':
    split3dfile=sys.argv[1]
    splitMEfile=sys.argv[2]
    splitMLfile=sys.argv[3]
    refsplit3d=sys.argv[4]
    refsplitME=sys.argv[5]
    refsplitML=sys.argv[6]
    outfile=sys.argv[7]

    IMD=import_split(split3dfile)
    ME=import_split(splitMEfile)
    ML=import_split(splitMLfile)
    rIMD=import_split(refsplit3d)
    rME=import_split(refsplitME)
    rML=import_split(refsplitML)

    ref_branches=comparing_references(rIMD, rME, rML)

    comparing_3d_ME_ML(IMD, ME, ML, ref_branches, outfile)
