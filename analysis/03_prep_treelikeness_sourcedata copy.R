library(phangorn)
library(ggplot2)
library(phangorn)

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
aligners = aligners[1]
trimming = trimming[1]

for (al in aligners){
    for (tr in trimming){

        df_complete = data.frame()
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

            #processing 3d dm
            rownames(dm3d)=dm3d$V1
            nm2=dm3d$V1
            dm3d=dm3d[,-1]
            colnames(dm3d)=nm2
            M3d=data.matrix(dm3d)
            m3d=c(M3d)

            #processing tm dm
            rownames(dmtm)=dmtm$V1
            nm3=dmtm$V1
            dmtm=dmtm[,-1]
            colnames(dmtm)=nm3
            Mtm=data.matrix(dmtm)
            mtm=c(Mtm)

            #computing ME patristic distance matrix
            P1me=cophenetic(tr1d)
            P1ml=cophenetic(trml)
            P3me=cophenetic(tr3d)
            P3tm=cophenetic(trtm)
            
            seqs =c()
            p1dme=c()
            d1dme=c()
            p1dml=c()
            p3dme=c()
            d3dme=c()
            p3tm=c()
            d3tm=c()

            for (i in seq_along(nm)){
                v_i = nm[i]
                for (j in seq_along(nm)){
                    v_j = nm[j]
                    if(i <= j){
                        next
                    }
                    seq   = paste(v_i, v_j, sep = "_")
                    seqs  = c(seqs, seq)
                    p1dme = c(p1dme, P1me[v_i,v_j])
                    p1dml = c(p1dml, P1ml[v_i,v_j])
                    d1dme = c(d1dme, M1d[v_i, v_j])
                    d3dme = c(d3dme, M3d[v_i, v_j])
                    p3dme = c(p3dme, P3me[v_i,v_j])
                    p3tm  = c(p3tm, P3tm[v_i,v_j])
                    d3tm  = c(d3tm, Mtm[v_i, v_j])
                }
            }

            # add in data frame
            df = data.frame(seqs, d1dme, p1dme, p1dml, d3dme, p3dme, p3tm, d3tm)
            colnames(df) = c('seqs', 'dME', 'pME', 'pML', 'dIMD', 'pIMD', 'pTM', 'dTM')
            # add family column
            df$family = fam
            df_complete = rbind(df_complete, df)
        }

        write.table(df_complete, file=paste(source_data,al,'_',tr,'_input_patristic_distances_raw.txt', sep=''), append=F,quote=F, sep = " ", eol = "\n", na = "NA", dec = ".",row.names =F, col.names =T)


    }
}

