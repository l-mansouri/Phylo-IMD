library(phangorn)
library(ggplot2)
library(ggExtra)

source_data = '/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data/'
fl=read.table('list_of_families_with_all_rep_in_3d')[,1]
setwd('/home/luisasantus/Desktop/crg_cluster/NF_draft/')


# -----------------------------------------------------------------------------
# Overview of the script: 
# read in trees and calculate shared branches
# -----------------------------------------------------------------------------

aligners=c('sap_tmalign', 'tcoffee')
trimming=c('untrimmed', 'trimmed')
        
for (al in aligners){
    print(al)
    for (tr in trimming){
        print(tr)
        TTB=c()
        TB=c()
        OKB_ML=c()
        OKB_ME=c()
        for (i in 1:100){
            total_total=c()
            total_branches=c()
            ok_branchesME=c()
            ok_branchesML=c()

            for (fam in fl){
                print(fam)
                #importing ML tree
                trml = read.tree(paste(al, '_ML_', tr,'_trees/',fam, '_', al, '_ML_', tr, '.nwk',sep=''))
                #importing 3d tree
                tr3d = read.tree(paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_',tr,'.nwk', sep=''))
                #tr3d = read.tree(paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_ME_',tr,'.nwk', sep=''))
                #importing ME
                tr1d = read.tree(paste(al, '_ME_', tr,'_trees/',fam, '_', al, '_ME_', tr, '.nwk',sep=''))

                b3d=tr3d$node.label
                b3d=b3d[-1]
                b3d=as.numeric(b3d)
                b3d[is.na(b3d)]<-0

                isinml=prop.clades(tr3d, trml)
                isinml[is.na(isinml)]=0
                isinml=isinml[-1]

                ttb=length(b3d)-1
                tb=length(b3d[b3d>=i])

                isinml=b3d*isinml
                okbml=length(isinml[isinml>=i])

                isinme=prop.clades(tr3d, tr1d)
                isinme[is.na(isinme)]=0
                isinme=isinme[-1]

                isinme=b3d*isinme
                okbme=length(isinme[isinme>=i])

                total_total=c(total_total,ttb)
                total_branches=c(total_branches, tb)
                ok_branchesML=c(ok_branchesML, okbml)
                ok_branchesME=c(ok_branchesME, okbme)
            }
            TTB=c(TTB,sum(total_total))
            TB=c(TB, sum(total_branches))
            OKB_ML=c(OKB_ML,sum(ok_branchesML))
            OKB_ME=c(OKB_ME,sum(ok_branchesME))
        }
        DF=data.frame(pos=1:100, OKB_ML, OKB_ME, TB, TTB)
        write.table(DF, file=paste(source_data,al,'_',tr,'counts_of_found_in_ML_and_ME_with_moving_IMD_BS_thr.txt', sep=''), row.names=F, col.names=T, quote=F)
    }
}



# -----------------------------------------------------------------------------
# Overview of the script: 
# read in trees and calculate RF distance
# -----------------------------------------------------------------------------


al="mTMalign"
tr="untrimmed"

avg_bs_ME=c()
avg_bs_3d=c()
avg_bs_ML=c()
RF_3d_ME=c()
RF_3d_ML=c()
for (fam in fl){
    #importing ME tree
    tr1d = read.tree(paste(al, '_ME_', tr,'_trees/',fam, '_', al, '_ME_', tr, '.nwk',sep=''))
    #importing ML tree
    trml = read.tree(paste(al, '_ML_', tr,'_trees/',fam, '_', al, '_ML_', tr, '.nwk',sep=''))
    #importing 3d tree
    #tr3d = read.tree(paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_',tr,'.nwk', sep=''))
    tr3d = read.tree(paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_ME_',tr,'.nwk', sep=''))
    #processing ME distance matrix
    b1d=tr1d$node.label
    b1d=b1d[-1]
    b1d=as.numeric(b1d)
    b1d[is.na(b1d)]<-0
    av1d=mean(b1d)

    b3d=tr3d$node.label
    b3d=b3d[-1]
    b3d=as.numeric(b3d)
    b3d[is.na(b3d)]<-0
    av3d=mean(b3d)

    bml=trml$node.label
    bml=bml[-1]
    bml=as.numeric(bml)
    bml[is.na(bml)]<-0
    avml=mean(bml)

    rf3dme=RF.dist(tr1d, tr3d, normalize=T)
    rf3dml=RF.dist(trml, tr3d, normalize=T)
    
    avg_bs_ME=c(avg_bs_ME, av1d)
    avg_bs_3d=c(avg_bs_3d, av3d )
    avg_bs_ML=c(avg_bs_ML, avml)
    RF_3d_ME=c(RF_3d_ME, rf3dme)
    RF_3d_ML=c(RF_3d_ML, rf3dml)
}
df=data.frame(
    family=fl,
    average_1d=as.numeric(avg_bs_ME),
    average_3d=as.numeric(avg_bs_3d),
    average_ML=as.numeric(avg_bs_ML),
    RF_3d_ME,
    RF_3d_ML
)
C3dme=round(cor(avg_bs_ME, avg_bs_3d),2)
C3dml=round(cor(avg_bs_ML, avg_bs_3d),2)
CrfME=round(cor(avg_bs_3d,RF_3d_ME ),2)
CrfML=round(cor(avg_bs_3d,RF_3d_ML ),2)


write.table(df, file = paste(source_data,"bootstrap_table.csv", sep = "/"), sep = ",", qmethod = "double", row.names = FALSE)
