
# ! PIPELINE OUTPUT FOLDER, CHANGE ACCORDINGLY
output_folder = '/home/luisasantus/Desktop/crg_cluster/NF_draft/'
#output_folder = "./results/"

# -----------------------------------------------------------------------------
# Overview of the script: 
# read in irmsd files and extract the total NiRMSD line
# save the NiRMSD in a df ready to be plotted
# -----------------------------------------------------------------------------


tcoffee_nirmsd <- list.files(path= paste(output_folder, "tcoffee_NiRMSD", sep = "/"), pattern = '*.NiRMSD', recursive = F, full.names = TRUE)
mtmalign_nirmsd <- list.files(path= paste(output_folder, "mTMalign_NiRMSD", sep = "/"), pattern = '*.NiRMSD', recursive = F, full.names = TRUE)
sap_tmalign_nirmsd <- list.files(path= paste(output_folder, "sap_tmalign_NiRMSD", sep = "/"), pattern = '*.NiRMSD', recursive = F, full.names = TRUE)
fl=read.table('source_data/list_of_families_with_all_rep_in_3d')[,1]

parse_nirmsd <- function(nirmsd_files, tag = "method"){
    # read in each file and then extract the total NiRMSD line and then save it in df 
    nirmsd_df <- data.frame()
    for (file in nirmsd_files){
        nirmsd <- readLines(file)
        nirmsd <- nirmsd[grep("TOTAL", nirmsd)]
        nirmsd <- nirmsd[grep("NiRMSD", nirmsd)]
        nirmsd <- gsub("TOTAL", "", nirmsd)
        nirmsd <- gsub("NiRMSD:", "", nirmsd)
        nirmsd <- gsub("Angs", "", nirmsd)
        nirmsd <- gsub(" ", "", nirmsd)
        nirmsd <- gsub("\t", "", nirmsd)
        nirmsd <- as.numeric(nirmsd)
        # extract name of the family and plsit _ and take the first element
        family <- strsplit(basename(file), "_")[[1]][1]
        # add nirmsd and fam to the df
        nirmsd_df <- rbind(nirmsd_df, data.frame(family, nirmsd))
    }
    nirmsd_df$type <- tag
    return (nirmsd_df)
}

table_tcoffee = parse_nirmsd(tcoffee_nirmsd, "T-Coffee")
table_mtmalign = parse_nirmsd(mtmalign_nirmsd, "mTM-align")
table_sap_tmalign = parse_nirmsd(sap_tmalign_nirmsd, "3D-Coffee")

df_tc=table_tcoffee[table_tcoffee$family %in% fl,]
df_mt=table_mtmalign[table_mtmalign$family %in% fl,]
df_st=table_sap_tmalign[table_sap_tmalign$family %in% fl,]

DF2=data.frame(NiRMSD=c(df_tc$nirmsd, df_st$nirmsd, df_mt$nirmsd),
              type=factor(rep(c('T-Coffee', '3D-Coffee', 'mTM-align'), each=(length(df_tc[,1]))), level= c('T-Coffee', '3D-Coffee', 'mTM-align'))
)
write.table(DF2, file = "source_data/S1_NiRMSD_table.csv", sep = ",", qmethod = "double", row.names = FALSE)
