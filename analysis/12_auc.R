# import file in same directory called utils_auc.R
# move wd to location of file
setwd("/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis")
source("utils_auc.R")
library(pROC)


#!/usr/bin/env Rscript
source_data = "/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data"
# OUTPUT DIRECTORY OF SPLIT FILES 
output_dir = "/home/luisasantus/Desktop/crg_cluster/newphylo/NF_TMalign/mTMalign/titration_every_5_bootstrap_200_columns/split_files"


# FAMILIES 
families =read.table(paste(source_data,'list_of_families_for_titration', sep = "/"))[,1]
# TREE TYPE TO BE EVALUATED 
tree_types = c("ME","ML", "IMD")
#tree_types = c("ME")

# REFERENCE DEFINITION 
reference = c("IMD", "ME", "IMD+ME", "ML", "IMD+ML", "ME+ML", "IMD+ME+ML")
# BOOTSTRAP THRESHOLD
bs_thresholds = c(0,80,100)
#bs_thresholds = c(80)

# EVALUATED BOOTSTRAP 
bs_combo = c("IMD","ME","ML", "ME+IMD", "ML+IMD", "ME+ML", "ME+ML+IMD")
# NUMBER OF COLUMNS
ncols = c(25,100,200)
#ncols = c(25)
# COMBINATION MODE 
comb_mode = c("arithmetic_average","geometric_average", "min","max")
#comb_mode = c("arithmetic_average")
auc_complete_df <- data.frame()
row_index <- 1
for (fam in families){
    print(fam)
    for (tree_type in tree_types){
        for (ref in reference){
            for (bs_threshold in bs_thresholds){
                # --------------------------------------------------
                #           EXTRACT PROVEN POSITIVES
                # --------------------------------------------------
                proven_positives = get_proven_positives(fam, tree_type, ref, bs_threshold)
                if(sum(proven_positives) == 0 | sum(proven_positives) == length(proven_positives)){
                    next
                }
                for (bs_type in bs_combo){
                    for (ncol in ncols){
                        # --------------------------------------------------
                        #           EXTRACT PREDICTED VALUES (BS)
                        # --------------------------------------------------

                        for (mode in comb_mode){
                            predicted_values = get_splits_bs(fam, tree_type, bs_type, ncol,mode)
                            # --------------------------------------------------
                            #                   COMPUTE AUC
                            # --------------------------------------------------                            
                            thresholds = seq(0, 100, 1)

                            auc_infos  = calculate_auc(proven_positives, predicted_values)
                            auc_v = auc_infos$auc
                            
                            # # calculate best threshold according to MCC
                            mcc_infos = get_best_mcc(predicted_values, proven_positives, thresholds)
                            best_threshold = mcc_infos$best_threshold
                            best_mcc = mcc_infos$best_mcc
                            sensitivity = mcc_infos$best_sensitivity
                            specificity = mcc_infos$best_specificity


                            # # add line to dataframe with all the infos
                            auc_complete_df <- rbind(auc_complete_df, data.frame(fam = fam, tree_type = tree_type, ref = ref, bs_threshold = bs_threshold, bs_type = bs_type, ncol = ncol, mode = mode, best_threshold = best_threshold, sensitivity = sensitivity, specificity = specificity, auc_v = auc_v, mcc = best_mcc))
                            row_index <- row_index + 1 
                        }
                    }
                }
            }
        }
    }
}


write.table(auc_complete_df, file = paste(source_data, "auc_complete_df_with_mcc.tsv", sep = "/"), sep = ",", row.names = FALSE)




