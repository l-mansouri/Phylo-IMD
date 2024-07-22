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
fl = fl[3]
setwd(output_folder)
aligners=c('mTMalign', 'sap_tmalign', 'tcoffee')
trimming=c('untrimmed', 'trimmed')

aligners = aligners[1]
trimming = trimming[1]

for (al in aligners){
    print(al)
    for (tr in trimming){

        COR1dme=c()
        IN1dme=c()
        PAT1dme=c()
        CORmeml=c()

        PAT1dml=c()
        PAT3dml=c()

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

            # -----------------------------------------------------------------------------
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
            # -----------------------------------------------------------------------------

            #processing ME distance matrix
            rownames(dm1d)=dm1d$V1
            nm=dm1d$V1
            dm1d=dm1d[,-1]
            colnames(dm1d)=nm
            M1d=data.matrix(dm1d)
            m1d=c(M1d)

            #computing ME patristic distance matrix
            P1me=cophenetic(tr1d)
            P1ml=cophenetic(trml)
            p1dme=c()
            seqs = c()
            for (i in seq_along(nm)){
                v_i = nm[i]
                for (j in seq_along(nm)){
                    v_j = nm[j]
                    if(i <= j){
                        next
                    }
                    seq = paste(v_i, v_j, sep = "_")
                    seqs = c(seqs, seq)
                    p1dme=c(p1dme, P1me[v_i,v_j])
                    p1dml=c(p1dml, P1ml[v_i,v_j])
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
            for (i in seq_along(nm2)){
                v_i = nm2[i]
                for (j in seq_along(nm2)){
                    v_j = nm2[j]
                    if(i <= j){
                        next
                    }
                    p3dme=c(p3dme, P3me[v_i,v_j])
                }
            }
            #compute patristic of ML for 3d correlations
            P3ml=P1ml
            p3dml=c()
            for (i in seq_along(nm2)){
                v_i = nm2[i]
                for (j in seq_along(nm2)){
                    v_j = nm2[j]
                    p3dml=c(p3dml, P3ml[v_i,v_j])
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
                for (i in seq_along(nm3)){
                    v_i = nm3[i]
                    v_j = nm3[j]
                    for (j in seq_along(nm3)){
                        if(i <= j){
                            next
                        }
                        p3tm=c(p3tm, P3tm[v_i,v_j])
                    }
                }
                #compute patristic of ML for 3d correlations
                Ptml=P1ml
                ptmml=c()
                for (i in seq_along(nm3)){
                    v_i = nm3[i]
                    for (j in seq_along(nm3)){
                        v_j = nm3[j]
                        ptmml=c(ptmml, Ptml[v_i,v_j])
                    }
                }
            }
            #collecting correlation, input, patristic, correlation to ML patristic
        #     COR1dme=c(COR1dme, cor(m1d, p1dme))
              IN1dme=c(IN1dme, m1d)
        #     PAT1dme=c(PAT1dme, p1dme)
        #     CORmeml=c(CORmeml, cor(m1d, p1dml))

        #     PAT1dml=c(PAT1dml, p1dml)
        #     PAT3dml=c(PAT3dml, p3dml)

        #     COR3dme=c(COR3dme, cor(m3d, p3dme))
        #     IN3dme=c(IN3dme, m3d)
        #     PAT3dme=c(PAT3dme, p3dme)
        #     COR3dml=c(COR3dml, cor(m3d, p3dml))
        #     if(al =='mTMalign' && tr == 'untrimmed'){
        #         CORtmme=c(CORtmme, cor(mtm, p3tm))
        #         INtmme=c(INtmme, mtm)
        #         PATtmme=c(PATtmme, p3tm)
        #         CORtmml=c(CORtmml, cor(mtm, ptmml))
        #     }
        }
        # correlations=data.frame(COR1dme, CORmeml, COR3dme, COR3dml)

        # Create a data frame with the input distances values and the patristic distances values 

        if(al =='mTMalign' && tr == 'untrimmed'){
            correlations_TM=cbind(correlations, CORtmme, CORtmml)
            #write.table(correlations_TM, file=paste(source_data,al,'_ME_',tr,'_saturation_correlations.txt', sep=''), append=F,quote=F, sep = " ", eol = "\n", na = "NA", dec = ".",row.names =F, col.names =T)
            # add column names
            input_and_patristic_distances = data.frame(IN1dme, IN3dme, INtmme, PAT1dme, PAT3dme, PATtmme, PAT1dml, PAT3dml)
            colnames(input_and_patristic_distances) = c('IN1dme', 'IN3dme', 'INtmme', 'PAT1dme', 'PAT3dme', 'PATtmme', 'PAT1dml', 'PAT3dml')
            # write the data frame to a file
            write.table(input_and_patristic_distances, file=paste(source_data,al,'_',tr,'_input_patristic_distances_test.txt', sep=''), append=F,quote=F, sep = " ", eol = "\n", na = "NA", dec = ".",row.names =F, col.names =T)
        }else{
            #write.table(correlations, file=paste(source_data,al,'_',tr,'_saturation_correlations_self_ML.txt', sep=''), append=F,quote=F, sep = " ", eol = "\n", na = "NA", dec = ".",row.names =F, col.names =T)
        }

    }
}
