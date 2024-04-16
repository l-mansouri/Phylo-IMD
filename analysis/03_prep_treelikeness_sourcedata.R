library(phangorn)
library(ggplot2)

# ! PIPELINE OUTPUT FOLDER, CHANGE ACCORDINGLY
output_folder = '/home/luisasantus/Desktop/crg_cluster/newphylo/NF_draft/'
# ! SOURCE DATA FOLDER, CHANGE ACCORDINGLY
source_data = '/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data/'


# -----------------------------------------------------------------------------
# Overview of the script: 
# read in trees and extract distances with cophenetic function from phangorn
# read in distance matrices and extract distances
# compute correlation between input distances and patristic distances
# -----------------------------------------------------------------------------
fl = read.table(paste(source_data,'list_of_families_with_all_rep_in_3d', sep = "/"))[,1]
setwd(output_folder)
aligners=c('mTMalign', 'sap_tmalign', 'tcoffee')
trimming=c('untrimmed', 'trimmed')



for (al in aligners){
    print(al)
    for (tr in trimming){

        COR1dme=c()
        IN1dme=c()
        PAT1dme=c()
        CORmeml=c()

        COR3dme=c()
        IN3dme=c()
        PAT3dme=c()
        COR3dml=c()

        CORtmme=c()
        INtmme=c()
        PATtmme=c()
        CORtmml=c()
        for (fam in fl){
            print(fam)
            #importing ME distance matrix and tree
            dm1d = read.table(paste(al,'_ME_',tr,'_matrix/',fam,'_', al, '_',tr,'_ME.matrix', sep=''), skip=1)
            tr1d = read.tree(paste(al, '_ME_', tr,'_trees/',fam, '_', al, '_ME_', tr, '.nwk',sep=''))
            #importing ML tree
            trml = read.tree(paste(al, '_ML_', tr,'_trees/',fam, '_', al, '_ML_', tr, '.nwk',sep=''))
            #importing 3d distance matrix and tree
            dm3d = read.table(paste(al,'_3d_ME_',tr,'_matrices/',fam,'_', al, '_3d_',tr,'_matrix.txt', sep=''), skip=1)
            if (al=='mTMalign'){
                tr3d = read.tree(paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_ME_',tr,'.nwk', sep=''))
            } else{
                tr3d = read.tree(paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_',tr,'.nwk', sep=''))
            }
            #importing tm score matrix and tree if it's computed
            if(al =='mTMalign' && tr == 'untrimmed'){
                trtm = read.tree(paste('mTMalign_TM_ME_trees/',fam,'_mTMalign_TM_ME_untrimmed.nwk',sep=''))
                dmtm = read.table(paste('mTMalign_matrix/',fam,'_mTMalign.matrix_4_fastme', sep=''), skip=1)
            }
            #processing ME distance matrix
            rownames(dm1d)=dm1d$V1
            nm=dm1d$V1
            dm1d=dm1d[,-1]
            colnames(dm1d)=nm
            M1d=data.matrix(dm1d)
            m1d=c(M1d)
            #computing ME patristic distance matrix
            P1me=cophenetic(tr1d)
            p1dme=c()
            for (i in nm){
                for (j in nm){
                    p1dme=c(p1dme, P1me[i,j])
                }
            }
            #computing ML patristic distance matrix for 1d corr
            P1ml=cophenetic(trml)
            p1dml=c()
            for (i in nm){
                for (j in nm){
                    p1dml=c(p1dml, P1ml[i,j])
                }
            }
            #processing 3d dm
            rownames(dm3d)=dm3d$V1
            nm2=dm3d$V1
            dm3d=dm3d[,-1]
            colnames(dm3d)=nm2
            M3d=data.matrix(dm3d)
            m3d=c(M3d)
            #computing 3d patristic distance matrix
            P3me=cophenetic(tr3d)
            p3dme=c()
            for (i in nm2){
                for (j in nm2){
                    p3dme=c(p3dme, P3me[i,j])
                }
            }
            #compute patristic of ML for 3d correlations
            P3ml=P1ml
            p3dml=c()
            for (i in nm2){
                for (j in nm2){
                    p3dml=c(p3dml, P3ml[i,j])
                }
            }
            if(al =='mTMalign' && tr == 'untrimmed'){
                rownames(dmtm)=dmtm$V1
                nm3=dmtm$V1
                dmtm=dmtm[,-1]
                colnames(dmtm)=nm3
                Mtm=data.matrix(dmtm)
                mtm=c(Mtm)
                #computing 3d patristic distance matrix
                P3tm=cophenetic(trtm)
                p3tm=c()
                for (i in nm3){
                    for (j in nm3){
                        p3tm=c(p3tm, P3tm[i,j])
                    }
                }
                #compute patristic of ML for 3d correlations
                Ptml=P1ml
                ptmml=c()
                for (i in nm3){
                    for (j in nm3){
                        ptmml=c(ptmml, Ptml[i,j])
                    }
                }
            }
            #collecting correlation, input, patristic, correlation to ML patristic
            COR1dme=c(COR1dme, cor(m1d, p1dme))
            IN1dme=c(IN1dme, m1d)
            PAT1dme=c(PAT1dme, p1dme)
            CORmeml=c(CORmeml, cor(m1d, p1dml))

            COR3dme=c(COR3dme, cor(m3d, p3dme))
            IN3dme=c(IN3dme, m3d)
            PAT3dme=c(PAT3dme, p3dme)
            COR3dml=c(COR3dml, cor(m3d, p3dml))
            if(al =='mTMalign' && tr == 'untrimmed'){
                CORtmme=c(CORtmme, cor(mtm, p3tm))
                INtmme=c(INtmme, mtm)
                PATtmme=c(PATtmme, p3tm)
                CORtmml=c(CORtmml, cor(mtm, ptmml))
            }
        }
        correlations=data.frame(COR1dme, CORmeml, COR3dme, COR3dml)
        if(al =='mTMalign' && tr == 'untrimmed'){
            correlations_TM=cbind(correlations, CORtmme, CORtmml)
            write.table(correlations_TM, file=paste(source_data,al,'_ME_',tr,'_saturation_correlations.txt', sep=''), append=F,quote=F, sep = " ", eol = "\n", na = "NA", dec = ".",row.names =F, col.names =T)
        }else{
            write.table(correlations, file=paste(source_data,al,'_',tr,'_saturation_correlations_self_ML.txt', sep=''), append=F,quote=F, sep = " ", eol = "\n", na = "NA", dec = ".",row.names =F, col.names =T)
        }

    }
}
