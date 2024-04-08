# import file in same directory called utils_auc.R
# move wd to location of file
source("utils_auc.R")
library(pROC)
library(dplyr)
library(tidyr)

#!/usr/bin/env Rscript
source_data = "/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data"
# OUTPUT DIRECTORY OF SPLIT FILES 
output_dir = "/home/luisasantus/Desktop/crg_cluster/newphylo/NF_TMalign/mTMalign/titration_every_5_bootstrap_200_columns/split_files"

# Read in the big AUC table
path_auc_table = paste(source_data, "auc_complete_df_with_mcc_me.tsv", sep = "/")
auc_complete_df = read.table(path_auc_table, header = TRUE, sep = ",", stringsAsFactors = FALSE)



# --------------------------------------------------
#                   TABLE 3 
# --------------------------------------------------

tree_type = "ML"
mode = "arithmetic_average"
bs_threshold = 80
ncol = 25

t3_data <- auc_complete_df[which(auc_complete_df$tree_type == tree_type & auc_complete_df$mode == mode & auc_complete_df$bs_threshold == bs_threshold & auc_complete_df$ncol == ncol),]
# round by 3 decimal places
t3_data$auc_v <- round(t3_data$auc_v, 3)
t3_data_grouped <- t3_data %>% group_by(ref, bs_type) %>% summarise(mean_auc = mean(auc_v))
col_order_bs = c("IMD", "ME", "ME+IMD", "ML", "ML+IMD", "ME+ML", "ME+ML+IMD")
col_order_ref = c("IMD", "ME", "IMD+ME", "ML", "IMD+ML", "ME+ML", "IMD+ME+ML")

# round the mean_auc to 3 decimal places
t3_data_grouped$mean_auc <- round(t3_data_grouped$mean_auc, 3)
# order the columns ref and bs_type by col_order
t3_data_grouped$ref <- factor(t3_data_grouped$ref, levels = col_order_ref)
t3_data_grouped$bs_type <- factor(t3_data_grouped$bs_type, levels = col_order_bs)
t3_matrix <- t3_data_grouped %>% spread(bs_type, mean_auc)
# write the table to a file
write.table(t3_matrix, file = paste(source_data, "../tables/Table3.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)



# --------------------------------------------------
#          SUPPLEMENTARY TABLE 4
# --------------------------------------------------
# create table with sd values ( add +/- and sd to each cell ) of the t3_data_summary
t3_data_grouped_sd <- t3_data %>% group_by(ref, bs_type) %>% summarise(std_auc = sd(auc_v))
t3_data_grouped_sd$std_auc <- round(t3_data_grouped_sd$std_auc, 2)
t3_data_grouped_sd$ref <- factor(t3_data_grouped_sd$ref, levels = col_order_ref)
t3_data_grouped_sd$bs_type <- factor(t3_data_grouped_sd$bs_type, levels = col_order_bs)
t3_matrix_sd <- t3_data_grouped_sd %>% spread(bs_type, std_auc)
write.table(t3_matrix_sd, file = paste(source_data, "../tables/S4.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)


# --------------------------------------------------
#          SUPPLEMENTARY TABLE 5
# --------------------------------------------------


# merge imd_col and me_col by fam, ref, bs_threshold, ncol, mode
do_wilcoxon_on_all_families <- function(t3_data, col1_id, col2_id){
    col1 = t3_data[t3_data$bs_type == col1_id,]
    col2 = t3_data[t3_data$bs_type == col2_id,]
    merged = merge(col1, col2, by = c("fam", "ref", "bs_threshold", "ncol", "mode"), all = TRUE)
    pval = wilcox.test(merged$auc_v.x, merged$auc_v.y, paired = TRUE, alternative = "greater")$p.value
    return(pval)
}


do_wilcoxon_on_merged_families <- function(t3_data_grouped, col1_id, col2_id){
    # order the columns ref and bs_type by col_order
    col1 = t3_data_grouped[t3_data_grouped$bs_type == col1_id,]$mean_auc
    col2 = t3_data_grouped[t3_data_grouped$bs_type == col2_id,]$mean_auc
    pval<- wilcox.test(col1, col2, paired = TRUE, alternative = "greater")$p.value
    return(pval)
}


# extract all possible pairs 
# get colnames from t3_data_grouped
colnames = unique(t3_data_grouped$bs_type)
# make a wilcox test for each pair of colnames 
df_wilcox = data.frame()
for (i in 1:length(colnames)){
    for (j in 1:length(colnames)){
        #print(colnames[i])
        #print(colnames[j])
        pval = do_wilcoxon_on_all_families(t3_data, colnames[i], colnames[j])
        pval_merged = do_wilcoxon_on_merged_families(t3_data_grouped, colnames[i], colnames[j])
        current_df = data.frame(firs = colnames[i], second = colnames[j], pval = pval, pval_merged = pval_merged)
        df_wilcox = rbind(df_wilcox, current_df) 
    }
}



# make it into a matrix
# df_wilcox_matrix = df_wilcox[,c("firs", "second", "pval_merged")] %>% spread(second, pval_merged)
# df_wilcox_matrix <- format(df_wilcox_matrix, scientific = TRUE)
# write.table(df_wilcox_matrix, file = paste(source_data, "../tables/S5_greater_merged.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)

df_wilcox_matrix = df_wilcox[,c("firs", "second", "pval")] %>% spread(firs, pval)
df_wilcox_matrix <- format(df_wilcox_matrix, scientific = TRUE)
df_wilcox_matrix
# write the table to a file
write.table(df_wilcox_matrix, file = paste(source_data, "../tables/S5_greater.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)



# order the columns ref and bs_type by col_order

# --------------------------------------------------
#          SUPPLEMENTARY TABLE 6
# --------------------------------------------------
ncol = 25
mode = "arithmetic_average"
bs_threshold = 80
tree_type = "ME"

table_4 <- auc_complete_df[which(auc_complete_df$tree_type == tree_type & auc_complete_df$mode == mode & auc_complete_df$bs_threshold == bs_threshold & auc_complete_df$ncol == ncol),]
table_4_grouped <- table_4 %>% group_by(ref, bs_type) %>% summarise(mean_auc = mean(auc_v))
col_order_bs = c("IMD", "ME", "ME+IMD", "ML", "ML+IMD", "ME+ML", "ME+ML+IMD")
col_order_ref = c("IMD", "ME", "IMD+ME", "ML", "IMD+ML", "ME+ML", "IMD+ME+ML")
table_4_grouped$mean_auc <- round(table_4_grouped$mean_auc, 3)
table_4_grouped$ref <- factor(table_4_grouped$ref, levels = col_order_ref)
table_4_grouped$bs_type <- factor(table_4_grouped$bs_type, levels = col_order_bs)
table_4_matrix <- table_4_grouped %>% spread(bs_type, mean_auc)
write.table(table_4_matrix, file = paste(source_data, "../tables/S6.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)

# do same with sd values
table_4_grouped_sd <- table_4 %>% group_by(ref, bs_type) %>% summarise(std_auc = sd(auc_v))
table_4_grouped_sd$std_auc <- round(table_4_grouped_sd$std_auc, 2)
table_4_grouped_sd$ref <- factor(table_4_grouped_sd$ref, levels = col_order_ref)
table_4_grouped_sd$bs_type <- factor(table_4_grouped_sd$bs_type, levels = col_order_bs)
table_4_matrix_sd <- table_4_grouped_sd %>% spread(bs_type, std_auc)
write.table(table_4_matrix_sd, file = paste(source_data, "../tables/S6_sd.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)



# --------------------------------------------------
#          SUPPLEMENTARY TABLE 7
# --------------------------------------------------
ncol = 25
tree_type = "ML"
bs_threshold = 80
ref  = "ME+ML"

s7 <- auc_complete_df[which(auc_complete_df$tree_type == tree_type & auc_complete_df$bs_threshold == bs_threshold & auc_complete_df$ncol == ncol & auc_complete_df$ref == ref),]
s7_grouped <- s7 %>% group_by(mode, bs_type) %>% summarise(mean_auc = mean(auc_v))
col_order = c("arithmetic_average", "geometric_average", "min", "max")
s7_grouped$mean_auc <- round(s7_grouped$mean_auc, 3)
s7_grouped$mode <- factor(s7_grouped$mode, levels = col_order)
col_order_bstype <- c("ME+IMD", "ML+IMD", "ME+ML", "ME+ML+IMD")
# filter out the bs_types that are not in col_order_bstype
s7_grouped <- s7_grouped[which(s7_grouped$bs_type %in% col_order_bstype),]
s7_matrix <- s7_grouped %>% spread(bs_type, mean_auc)
write.table(s7_matrix, file = paste(source_data, "../tables/S7.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)

# --------------------------------------------------
#          SUPPLEMENTARY TABLE 8
# --------------------------------------------------

# with sensitivity
t3_data_grouped_sensitivity <- t3_data %>% group_by(ref, bs_type) %>% summarise(mean_sensitivity = mean(sensitivity))
t3_data_grouped_sensitivity$mean_sensitivity <- round(t3_data_grouped_sensitivity$mean_sensitivity, 2)
t3_data_grouped_sensitivity$ref <- factor(t3_data_grouped_sensitivity$ref, levels = col_order_ref)
t3_data_grouped_sensitivity$bs_type <- factor(t3_data_grouped_sensitivity$bs_type, levels = col_order_bs)
s8_matrix <- t3_data_grouped_sensitivity %>% spread(bs_type, mean_sensitivity)
write.table(s8_matrix, file = paste(source_data, "../tables/S8_sens.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)

# with specificity
t3_data_grouped_specificity <- t3_data %>% group_by(ref, bs_type) %>% summarise(mean_specificity = mean(specificity))
t3_data_grouped_specificity$mean_specificity <- round(t3_data_grouped_specificity$mean_specificity, 2)
t3_data_grouped_specificity$ref <- factor(t3_data_grouped_specificity$ref, levels = col_order_ref)
t3_data_grouped_specificity$bs_type <- factor(t3_data_grouped_specificity$bs_type, levels = col_order_bs)
s8_matrix <- t3_data_grouped_specificity %>% spread(bs_type, mean_specificity)
write.table(s8_matrix, file = paste(source_data, "../tables/S8_spec.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)




# --------------------------------------------------
#          SUPPLEMENTARY TABLE 9
# --------------------------------------------------
t3_data_grouped_threshold <- t3_data %>% group_by(ref, bs_type) %>% summarise(best_threshold = mean(best_threshold))
t3_data_grouped_threshold$best_threshold <- round(t3_data_grouped_threshold$best_threshold, 0)
t3_data_grouped_threshold$ref <- factor(t3_data_grouped_threshold$ref, levels = col_order_ref)
t3_data_grouped_threshold$bs_type <- factor(t3_data_grouped_threshold$bs_type, levels = col_order_bs)
s9_matrix <- t3_data_grouped_threshold %>% spread(bs_type, best_threshold)
write.table(s9_matrix, file = paste(source_data, "../tables/S9_mean.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)

t3_data_grouped_threshold_sd <- t3_data %>% group_by(ref, bs_type) %>% summarise(std_threshold = sd(best_threshold))
t3_data_grouped_threshold_sd$std_threshold <- round(t3_data_grouped_threshold_sd$std_threshold, 0)
t3_data_grouped_threshold_sd$ref <- factor(t3_data_grouped_threshold_sd$ref, levels = col_order_ref)
t3_data_grouped_threshold_sd$bs_type <- factor(t3_data_grouped_threshold_sd$bs_type, levels = col_order_bs)
s9_matrix_sd <- t3_data_grouped_threshold_sd %>% spread(bs_type, std_threshold)
write.table(s9_matrix_sd, file = paste(source_data, "../tables/S9_sd.csv", sep = "/"), sep = ",", quote = FALSE, row.names = FALSE)
