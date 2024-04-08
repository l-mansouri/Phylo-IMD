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
tree_types = c("ML")

# REFERENCE DEFINITION 
reference = c("IMD","ME","ML", "IMD+ME", "IMD+ML", "ME+ML", "IMD+ME+ML")
reference = c("IMD")
# BOOTSTRAP THRESHOLD
bs_thresholds = c(0,80,100)
bs_thresholds = c(80)

# EVALUATED BOOTSTRAP 
bs_combo = c("IMD","ME","ML", "IMD+ME", "IMD+ML", "ME+ML", "IMD+ME+ML")
bs_combo = c("ME")

# NUMBER OF COLUMNS
ncols = c(25,100,200)
ncols = c(25)
# COMBINATION MODE 
comb_mode = c("arithmetic_average","geometric_average", "min","max")
comb_mode = c("arithmetic_average")

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

                            roc_data = roc(proven_positives, predicted_values,ci.thresholds = thresholds)
                            auc_value = as.numeric(auc(roc_data))

                            # # calculate best threshold according to MCC
                            best_threshold = 0
                            best_mcc = -1
                            auc_v = 0 
                            prev_tpr = 0
                            prev_fpr = 0
                            for( threshold in thresholds){
                                predicted_labels <- ifelse(predicted_values > threshold, 1, 0)
                                predicted_labels = factor(predicted_labels, levels = c(0,1))

                                confusion_matrix = table(predicted_labels, proven_positives)
                                tp = confusion_matrix[2,2]
                                tn = confusion_matrix[1,1]
                                fp = confusion_matrix[1,2]
                                fn = confusion_matrix[2,1]
                                tpr = tp / (tp + fn)
                                fpr = fp / (fp + tn)
                                # if they are none, set them to 0
                                if(is.nan(tpr)){
                                    tpr = 0
                                }
                                if(is.nan(fpr)){
                                    fpr = 0
                                }

                                mcc = (tp*tn - fp*fn) / sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
                                # calculate auc with mcc
                                if (fp == 0 & tp > 0){
                                    prev_tpr = tpr
                                }
                                else if(fp > 0 & tp > 0){
                                    auc_v = auc_v + ((fpr + prev_fpr) * (tpr + prev_tpr) / 2)
                                }
                                
                                if(is.nan(mcc)){
                                    mcc = -1
                                }
                                if(mcc > best_mcc){
                                    best_mcc = mcc
                                    best_threshold = threshold
                                }
                            }
                            # # calculate sensitivity and specificity using the threshold best_threshold
                            predicted_labels = ifelse(predicted_values > best_threshold, 1, 0)
                            predicted_labels = factor(predicted_labels, levels = c(0,1))
                            confusion_matrix = table(predicted_labels, proven_positives)
                            sensitivity = confusion_matrix[2,2] / sum(proven_positives)
                            specificity = confusion_matrix[1,1] / sum(proven_positives == 0)

                            # # add line to dataframe with all the infos
                            auc_complete_df <- rbind(auc_complete_df, data.frame(fam = fam, tree_type = tree_type, ref = ref, bs_threshold = bs_threshold, bs_type = bs_type, ncol = ncol, mode = mode, auc_value = auc_value, best_threshold = best_threshold, sensitivity = sensitivity, specificity = specificity, auc_v = auc_v, mcc = best_mcc))
                            row_index <- row_index + 1 
                        }
                    }
                }
            }
        }
    }
}


mean(auc_complete_df$auc_value)

mean(round(auc_complete_df$auc_value,2))

mean(auc_complete_df$auc_v)


df <- auc_complete_df

df$auc_v

print(mb)
mean(auc_complete_df$auc_value)

# Save the dataframe
write.table(auc_complete_df, file = paste(source_data, "auc_complete_df_with_mcc_small.tsv", sep = "/"), sep = ",", row.names = FALSE)


calculate_auc <- function(true_positives, predicted_positives) {
    # sort the predicted positives in descending order (and keep the true positives in the same order)
    n <- length(true_positives)
    print(n)
    data <- data.frame(score = predicted_positives, label = true_positives)
    print(data)
    data <- data[order(-data$score), ]

    auc <- 0.0
    prev_fpr <- 0.0
    prev_tpr <- 0.0
    tp <- 0
    fp <- 0
    positives <- 0
    negatives <- 0

    for (i in 1:n) {
        if (data[i,]$label == 1) {
            positives <- positives + 1
        } else {
            negatives <- negatives + 1
        }
    }

    i <- 1
    while (i <= n) {
    j <- i
    tie_fp <- 0
    tie_tp <- 0
    
    # Handle ties by aggregating them
    while (j <= n && data[j,]$score == data[i,]$score) {
        if (data[j,]$label == 1) {
        tie_tp <- tie_tp + 1
        } else {
        tie_fp <- tie_fp + 1
        }
        j <- j + 1
    }
    
    tp <- tp + tie_tp
    fp <- fp + tie_fp
    
    tpr <- tp / positives
    fpr <- fp / negatives
    
    if (fp == 0 && tp > 0) {
        prev_tpr <- tpr
    } else if (i > 1 || (fp > 0 && tp > 0)) {
        auc <- auc + (fpr - prev_fpr) * (tpr + prev_tpr) / 2.0
        prev_fpr <- fpr
        prev_tpr <- tpr
    }
    
    i <- j  # Move to the next group of scores
    }
  
  return(auc)
}

# Example usage
true_positives <- c(1, 0, 1, 1, 0)
predicted_positives <- c(0.9, 0.7, 0.6, 0.4, 0.2)

auc <- calculate_auc(true_positives, predicted_positives)
print(auc)
