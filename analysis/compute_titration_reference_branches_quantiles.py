#!/usr/bin/env python
import sys
import numpy as np
import os


def read_RF_for_fraction_2(frac1d1d,frac3d3d,filename):
    with open(filename) as f:
        for line in f:
            line=line.rstrip().split(' ')
            frac1d1d.setdefault(line[0],[]).append(float(line[1]))
            frac3d3d.setdefault(line[0],[]).append(float(line[2]))
    return frac1d1d,frac3d3d

def read_RF_for_fraction_3(frac1d1d,frac3d3d,fracmdmd, filename):
    with open(filename) as f:
        for line in f:
            line=line.rstrip().split(' ')
            frac1d1d.setdefault(line[0],[]).append(float(line[1]))
            frac3d3d.setdefault(line[0],[]).append(float(line[2]))
            fracmdmd.setdefault(line[0],[]).append(float(line[3]))
    return frac1d1d,frac3d3d, fracmdmd


if __name__ == '__main__':
    # read in arguments 
    lista=sys.argv[1]
    titration_dir=str(sys.argv[2])
    decile = lista.split('/')[-1].split('.txt')[0]
    # read in family list
    l=[]
    for line in open(lista):
        line=line.rstrip()
        l.append(line)
    outfilename=decile+'_titration_every_5_reference_branches_ME+3d+ML.txt'
    outfilename = os.path.join("source_data/", outfilename)
    frac1d1d={}
    frac3d3d={}
    fracmdmd={}
    cnt={}
    for family in l:
        filename=str(titration_dir)+'/ref_br_ME+3d+ML/'+str(family)+'_avg_fr_ref_br_per_fraction_shared_1d_3d_ML'
        frac1d1d,frac3d3d,fracmdmd=read_RF_for_fraction_3(frac1d1d,frac3d3d,fracmdmd,filename)
    out=open(outfilename, 'w')
    #ks=frac.keys()
    for k in frac1d1d.keys():
        vals1d1d=np.array(frac1d1d[k])
        avg1d1d=np.mean(vals1d1d)
        stdev1d1d=np.std(vals1d1d)
        vals3d3d=np.array(frac3d3d[k])
        avg3d3d=np.mean(vals3d3d)
        stdev3d3d=np.std(vals3d3d)
        valsmdmd=np.array(fracmdmd[k])
        avgmdmd=np.mean(valsmdmd)
        stdevmdmd=np.std(valsmdmd)
        out.write(str(k)+'\t'+str(avg1d1d)+'\t'+str(stdev1d1d)+'\t'+str(avg3d3d)+'\t'+str(stdev3d3d)+'\t'+str(avgmdmd)+'\t'+str(stdevmdmd)+'\t'+'100'+'\n')
    out.close()

