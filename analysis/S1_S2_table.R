library(phangorn)

aligners=c('mTMalign', "sap_tmalign", "tcoffee")
trimming=c('untrimmed', "trimmed")

# ! PIPELINE OUTPUT FOLDER, CHANGE ACCORDINGLY
output_folder = '/home/luisasantus/Desktop/crg_cluster/newphylo/NF_draft/'
# ! SOURCE DATA FOLDER, CHANGE ACCORDINGLY
source_data = '/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data/'



fl = read.table(paste(source_data,'/list_of_families_with_all_rep_in_3d', sep = "/"), header=F)[,1]
tables = paste( source_data, "../tables/", sep = "/" )
setwd(output_folder)

complete_df = data.frame()
for (al in aligners){    
    for (tr in trimming){

        EV_ME=c()
        EV_3d=c()
        bs_ME=c()
        bs_ML=c()
        bs_3d=c()
        RF_3d_ME=c()
        RF_3d_ML=c()
        RF_ME_ML=c()
        RF_3d_ME_75=c()
        RF_3d_ML_75=c()
        RF_ME_ML_75=c()

        SB_ME_ML=c()
        SB_3D_ML=c()
        SB_TM_ML=c()

        SB_ME_ML_100=c()
        SB_3D_ML_100=c()

        TOT_3D_100=c()
        TOT_ML_100=c()
        TOT_ME_100=c()
        TOTAL_BRANCHES=c()
        #additional vectors for the TM trees
        if(al =='mTMalign' && tr == 'untrimmed'){
            EV_TM=c()
            RF_TM_ME=c()
            RF_TM_ML=c()
        }
        for (fam in fl){
            print(fam)
            #importing ME tree dm
            dm1d = read.table(paste(al,'_ME_',tr,'_matrix/',fam,'_', al, '_',tr,'_ME.matrix', sep=''), skip=1)
            tr1d = read.tree(paste(al, '_ME_', tr,'_trees/',fam, '_', al, '_ME_', tr, '.nwk',sep=''))
            #importing ML tree
            trml = read.tree(paste(al, '_ML_', tr,'_trees/',fam, '_', al, '_ML_', tr, '.nwk',sep=''))
            #importing 3d tree dm
            dm3d = read.table(paste(al,'_3d_ME_',tr,'_matrices/',fam,'_', al, '_3d_',tr,'_matrix.txt', sep=''), skip=1)
            if (al=='mTMalign'){
                tr3d = read.tree(paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_ME_',tr,'.nwk', sep=''))
            } else{
                tr3d = read.tree(paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_',tr,'.nwk', sep=''))
            }
            #importing TM trees dm
            if(al =='mTMalign' && tr == 'untrimmed'){
                trtm = read.tree(paste('mTMalign_TM_ME_trees/',fam,'_mTMalign_TM_ME_untrimmed.nwk',sep=''))
                dmtm = read.table(paste('mTMalign_matrix/',fam,'_mTMalign.matrix_4_fastme', sep=''), skip=1)
            }
            #computing EV
            #processing ME distance matrix
            rownames(dm1d)=dm1d$V1
            nm=dm1d$V1
            dm1d=dm1d[,-1]
            colnames(dm1d)=nm
            M1d=data.matrix(dm1d)
            #computing ME patristic distance matrix
            P1me=cophenetic(tr1d)
            #processing 3d dm
            rownames(dm3d)=dm3d$V1
            nm2=dm3d$V1
            dm3d=dm3d[,-1]
            colnames(dm3d)=nm2
            M3d=data.matrix(dm3d)
            #computing 3d patristic distance matrix
            P3me=cophenetic(tr3d)

            #processing bootstrap
            b1d=tr1d$node.label
            b1d=b1d[-1]
            b1d=as.numeric(b1d)
            b1d[is.na(b1d)]<-0
            bs_ME=c(bs_ME, b1d)

            b3d=tr3d$node.label
            b3d=b3d[-1]
            b3d=as.numeric(b3d)
            b3d[is.na(b3d)]<-0
            bs_3d=c(bs_3d, b3d)

            bml=trml$node.label
            bml=bml[-1]
            bml=as.numeric(bml)
            bml[is.na(bml)]<-0
            bs_ML=c(bs_ML, bml)

            # select trees with avg_bs > 75
            avg_bs_ME = mean(b1d)
            avg_bs_3d = mean(b3d)
            avg_bs_ML = mean(bml)
            if (avg_bs_ME >= 75){
                rfmeml75 = RF.dist(tr1d, trml, normalize=T, rooted=F)
            }
            if (avg_bs_3d >= 75){
                rf3dml75 = RF.dist(tr3d, trml, normalize=T, rooted=F)
            }

            #shared branches
            

            sbmeml = prop.clades(tr1d, trml)
            sbmeml <- sbmeml[-1]
            sbmeml[is.na(sbmeml)]<-0
            SB_ME_ML=c(SB_ME_ML, sum(sbmeml))
            TOTAL_BRANCHES = c(TOTAL_BRANCHES, length(sbmeml))

            sb3dml = prop.clades(tr3d, trml)
            sb3dml <- sb3dml[-1]
            sb3dml[is.na(sb3dml)]<-0
            SB_3D_ML=c(SB_3D_ML, sum(sb3dml))

            # check if bs is == 100, return vector with 1 or 0
            b3d_100 = ifelse(b3d == 100, 1, 0)
            TOT_3D_100 = c(TOT_3D_100, sum(b3d_100))

            b1d_100 = ifelse(b1d == 100, 1, 0)
            TOT_ME_100 = c(TOT_ME_100, sum(b1d_100))
            # calculate how many branches are shared with ML tree
            sbmeml_100 = sbmeml * b1d_100
            SB_ME_ML_100 = c(SB_ME_ML_100, sum(sbmeml_100))

            sb3dml_100 = sb3dml * b3d_100
            SB_3D_ML_100 = c(SB_3D_ML_100, sum(sb3dml_100))

            #computing RFs
            rf3dme=RF.dist(tr3d, tr1d, normalize=T, rooted = F)
            rf3dml=RF.dist(tr3d, trml, normalize=T, rooted = F)
            rfmeml=RF.dist(tr1d, trml, normalize=T, rooted = F)

            #collecting RFs
            if(al =='mTMalign' && tr == 'untrimmed'){
                rftmme=RF.dist(trtm, tr1d, normalize=T, rooted=F)
                rftmml=RF.dist(trtm, trml, normalize=T, rooted=F)
                RF_TM_ME=c(RF_TM_ME, rftmme)
                RF_TM_ML=c(RF_TM_ML, rftmml)

                sbtmml = prop.clades(trtm, trml)
                sbtmml <- sbtmml[-1]
                sbtmml[is.na(sbtmml)]<-0
                SB_TM_ML=c(SB_TM_ML, sum(sbtmml))
            }

            RF_3d_ME=c(RF_3d_ME, rf3dme)
            RF_3d_ML=c(RF_3d_ML, rf3dml)
            RF_ME_ML=c(RF_ME_ML, rfmeml)
            if (avg_bs_ME >= 75){
                RF_ME_ML_75=c(RF_ME_ML_75, rfmeml75)
            }
            if (avg_bs_3d >= 75){
                RF_3d_ML_75=c(RF_3d_ML_75, rf3dml75)
            }

        }


        re_ME=c(al, tr, round(mean(bs_ME),2), round(mean(RF_ME_ML), 2), round(mean(RF_ME_ML_75), 2), round(mean(SB_ME_ML*100/TOTAL_BRANCHES),2), sum(SB_ME_ML), sum(SB_ME_ML_100), round(sum(SB_ME_ML_100)*100/sum(TOT_ME_100),2))
        re_3d=c(al, tr, round(mean(bs_3d),2), round(mean(RF_3d_ML),2), round(mean(RF_3d_ML_75),2), round(mean(SB_3D_ML*100/TOTAL_BRANCHES),2), sum(SB_3D_ML), sum(SB_3D_ML_100), round(sum(SB_3D_ML_100)*100/sum(TOT_3D_100),2))
        # Prep final dataframe
        df=as.data.frame(rbind(re_ME, re_3d))
        colnames(df)=c("aligner", "trimmer",'avg_BS','RFxML', 'RFxML_75', "Shared branches w/ ML", "Count shared branches", "Count 100% BS", "% 100% BS shared")
        
        if(al =='mTMalign' && tr == 'untrimmed'){
            re_TM=c(al, tr, NA, round(mean(RF_TM_ML),2), NA, round(mean(SB_TM_ML*100/TOTAL_BRANCHES),2), sum(SB_TM_ML), NA, NA,NA)
            df=rbind(df, re_TM)
            rownames(df)=c('Seq-ME','IMD-ME', 'TM-ME') 
        }else{
            rownames(df)=c('Seq-ME','IMD-ME')
        }
        complete_df = rbind(complete_df, df)
    }
}

# add congruence column (1-RF)
# make numeric
complete_df$RFxML <- as.numeric(as.character(complete_df$RFxML))
complete_df$congruence_ML <- 1-complete_df$RFxML
complete_df$RFxML_75 <- as.numeric(as.character(complete_df$RFxML_75))
complete_df$congruence_ML_75 <- 1-complete_df$RFxML_75

# ---------------------------------------
# Put together in the right format
# ---------------------------------------
complete_df
df <- complete_df
# select only the columns of interest
columns <- c("avg_BS","congruence_ML", "congruence_ML_75", "Shared branches w/ ML", "Count shared branches","Count 100% BS", "% 100% BS shared")
df_clean <- df[,columns]
# add trimmer and aligner as the first columns
df_clean$aligner <- df$aligner
df_clean$trimmer <- df$trimmer
df_clean$method <- rownames(df_clean)
df_clean$method <- gsub("[0-9]", "", df_clean$method)
df_clean <- df_clean[,c(8,9,10,1:7)]
# keep only seqme, tm me and imd me rows (in this order)
columns <- c("aligner", "trimming", "method","avg_BS","congruence_ML", "congruence_ML_75", "Shared branches w/ ML", "Count shared branches", "Count 100% BS", "% 100% BS shared")
df_clean = df_clean[c(1,3,2),columns]
df_clean
# rename columns
colnames(df_clean) <- c("aligner", "trimming", "method","Average BS","Congruence with Seq-ML trees", "Congruence with Seq-ML trees (BS>=75)", "Seq-ML Shared branches (%)", "Seq-ML Shared branches (count)", "BS=100 Seq-ML Shared branches (count)", "BS=100 Seq-ML Shared branches (%)")
write.table(df_clean, file=paste(tables, paste("Table1.csv", sep=''), sep = "/"), append=F,quote=F, sep = ",", eol = "\n", na = "NA", dec = ".",row.names =F, col.names =T)

# -----------------------------
# Extract patristic distances
# -----------------------------
setwd(source_data) 
mtm_untrimmed = read.table('mTMalign_ME_untrimmed_saturation_correlations.txt', header=T)
mtm_trimmed = read.table('mTMalign_trimmed_saturation_correlations_self_ML.txt', header=T)
sap_tm_untrimmed = read.table('sap_tmalign_untrimmed_saturation_correlations_self_ML.txt', header=T)
sap_tm_trimmed = read.table('sap_tmalign_trimmed_saturation_correlations_self_ML.txt', header=T)
tcoffee_untrimmed = read.table('tcoffee_untrimmed_saturation_correlations_self_ML.txt', header=T)
tcoffee_trimmed = read.table('tcoffee_trimmed_saturation_correlations_self_ML.txt', header=T)


calc_medians <- function(df, method, trimming){
    median_seq = round(median((df$COR1dme)**2),3)
    median_imd = round(median((df$COR3dme)**2),3)
    if((method == "mTMalign") & (trimming == "untrimmed")){
        median_tm  = round(median((df$CORtmme)**2),3)
    } else{
        median_tm = NA
    }
    df1 <- data.frame("aligner" = method, "trimming" = trimming, "method" = "Seq-ME", "R squared" = median_seq ) 
    df2 <- data.frame("aligner" = method, "trimming" = trimming, "method" = "IMD-ME", "R squared" = median_imd )
    df3 <- data.frame("aligner" = method, "trimming" = trimming, "method" = "TM-ME", "R squared" = median_tm )
    return(rbind(df1, df2, df3))
}


# dataframe with colnames , method, trimming, Seq-ME, IMD-ME, TM-ME
df_patristic = data.frame("aligner" = character(), "trimming" = character(), "Seq-ME" = numeric(), "IMD-ME" = numeric(), "TM-ME" = numeric())
df_patristic = rbind(df_patristic, calc_medians(mtm_untrimmed, "mTMalign", "untrimmed"))
df_patristic = rbind(df_patristic, calc_medians(mtm_trimmed, "mTMalign", "trimmed"))
df_patristic = rbind(df_patristic, calc_medians(sap_tm_untrimmed, "sap_tmalign", "untrimmed"))
df_patristic = rbind(df_patristic, calc_medians(sap_tm_trimmed, "sap_tmalign", "trimmed"))
df_patristic = rbind(df_patristic, calc_medians(tcoffee_untrimmed, "tcoffee", "untrimmed"))
df_patristic = rbind(df_patristic, calc_medians(tcoffee_trimmed, "tcoffee", "trimmed"))
# substitue . for - in colnames
colnames(df_patristic) <- gsub("\\.", "-", colnames(df_patristic))

df_patristic
# only untrimmed
df_patristic_untrimmed = df_patristic[df_patristic$trimming == "untrimmed",]
df_patristic_untrimmed

df_patristic_trimmed = df_patristic[df_patristic$trimming == "trimmed",]
df_patristic_trimmed
# -----------------------------
# Merge tables
# -----------------------------
df_clean_with_patristic = merge(df_clean, df_patristic, by=c("aligner", "trimming", "method"), all=T)
df_clean_with_patristic <- df_clean_with_patristic[,c(1,2,3,11,4:10)]
order_aligner = c("mTMalign", "sap_tmalign", "tcoffee")
order_trimming = c("untrimmed", "trimmed")
df_clean_with_patristic$aligner <- factor(df_clean_with_patristic$aligner, levels = order_aligner)
df_clean_with_patristic$trimming <- factor(df_clean_with_patristic$trimming, levels = order_trimming)
df_clean_with_patristic$method <- factor(df_clean_with_patristic$method, levels = c("Seq-ME","TM-ME","IMD-ME"))
df_clean_with_patristic

# -----------------------------------
# s1
# -----------------------------------

# sort by aligner, trimming and method with given order (order by aligner, then by trimming and then by method)
df_clean_with_patristic <- df_clean_with_patristic[order(df_clean_with_patristic$aligner, df_clean_with_patristic$trimming, df_clean_with_patristic$method),]
df_clean_with_patristic_untrimmed = df_clean_with_patristic[df_clean_with_patristic$trimming == "untrimmed",]
# save table
write.table(df_clean_with_patristic_untrimmed, file=paste(tables, paste("S1.csv", sep=''), sep = "/"), append=F,quote=F, sep = ",", eol = "\n", na = "NA", dec = ".",row.names =F, col.names =T)

# -----------------------------------
# s2
# -----------------------------------
# now trimmed
df_clean_with_patristic_trimmed = df_clean_with_patristic[df_clean_with_patristic$trimming == "trimmed",]
# write table
write.table(df_clean_with_patristic_trimmed, file=paste(tables, paste("S2.csv", sep=''), sep = "/"), append=F,quote=F, sep = ",", eol = "\n", na = "NA", dec = ".",row.names =F, col.names =T)


