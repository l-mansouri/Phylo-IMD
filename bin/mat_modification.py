import sys

def in_fasta(fastafile):
    ID=[]
    for l in open(fastafile):
        l=l.rstrip()
        if len(l)>0:
            if l[0]=='>':
                ID.append(l[1:])
    return ID

def in_mat(matrixfile,outfile, ID):
    out=open(outfile,'w')
    acc=0
    for l in open(matrixfile):
        l=l.rstrip()
        acc=acc+1
        if acc>1:
            out.write(str(ID[acc-2][0:10])+str(l[10:])+'\n')
        else:
            out.write(str(l)+'\n')
    out.close()

if __name__ == '__main__':
    fastafile=sys.argv[1]
    matrixfile=sys.argv[2]
    outfile=sys.argv[3]
# fastafile='mTMalign_matrix/PF00051_mTMalign.fa'
# matrixfile='mTMalign_matrix/PF00051_mTMalign.matrix'
# outfile='prova'
    ID=in_fasta(fastafile)
    in_mat(matrixfile,outfile, ID)
