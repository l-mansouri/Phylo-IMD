# import file in same directory called utils_auc.R
# move wd to location of file
setwd("/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis")
source("utils_auc.R")
library(pROC)
library(patchwork)

#!/usr/bin/env Rscript
source_data = "/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data"
# OUTPUT DIRECTORY OF SPLIT FILES 
output_dir = "/home/luisasantus/Desktop/crg_cluster/newphylo/NF_TMalign/mTMalign/titration_every_5_bootstrap_200_columns/split_files"

# FAMILIES 
families =read.table(paste(source_data,'list_of_families_for_titration', sep = "/"))[,1]


get_all_pp_and_pred <- function(tree_type, references, bs_threshold, bs_type, ncol, mode){
    proven_positives = c()
    predicted_values = c()
    for (fam in families){
        for (ref in references){
            print(ref)
            new_proven_positives = get_proven_positives(fam, tree_type, ref, bs_threshold)
            print("proved positives")
            if(sum(new_proven_positives) == 0 | sum(new_proven_positives) == length(new_proven_positives)){
                next
            }
            proven_positives = c(proven_positives, new_proven_positives)

            new_predicted_values = get_splits_bs(fam, tree_type, bs_type, ncol,mode)
            predicted_values = c(predicted_values, new_predicted_values)
        }

    }
    return(data.frame(proven_positives = proven_positives, predicted_values = predicted_values))
}

prep_df_roc <- function(tree_type, ref, bs_threshold, bs_type1, bs_type2, ncol, mode){
    df1 <- get_all_pp_and_pred(tree_type = tree_type, ref = ref, bs_threshold = bs_threshold, bs_type = bs_type1, ncol = ncol, mode = mode)
    df2 <- get_all_pp_and_pred(tree_type = tree_type, ref = ref, bs_threshold = bs_threshold, bs_type = bs_type2, ncol = ncol, mode = mode)
    roc1 <- roc(df1$proven_positives, df1$predicted_values)
    roc2 <- roc(df2$proven_positives, df2$predicted_values)

    roc1_df <- data.frame(roc1$thresholds, roc1$sensitivities, roc1$specificities)
    roc2_df <- data.frame(roc2$thresholds, roc2$sensitivities, roc2$specificities)
    colnames(roc1_df) <- c("threshold", "sensitivity", "specificity")
    colnames(roc2_df) <- c("threshold", "sensitivity", "specificity")
    roc1_df$roc <- bs_type1
    roc2_df$roc <- bs_type2
    roc_df <- rbind(roc1_df, roc2_df)
    roc_df$FPR <- 1 - roc_df$specificity

    # add auc values 
    auc1  <- auc(roc1)
    auc2  <- auc(roc2)
    roc_df$auc <- NA
    roc_df$auc[roc_df$roc == bs_type1] <- auc1
    roc_df$auc[roc_df$roc == bs_type2] <- auc2
    return(roc_df)

}

plot_ROC <- function(roc_df, path, title = "ROC curve" ){
    p1 <- ggplot(roc_df, aes(x=FPR, y=sensitivity, color=roc, linetype=roc)) +
        geom_line(lwd = 0.8) +
        geom_abline(intercept = 0, slope = 1, lty = 2) +
        theme(legend.position="bottomright") +
        ggtitle(title) +
        xlab("1 - specificity") +
        ylab("sensitivity") +
        scale_color_manual(values=c("blue", "red")) +
        scale_linetype_manual(values=c("solid", "solid")) +
        scale_x_continuous(limits = c(0,1)) +
        scale_y_continuous(limits = c(0,1)) +
        theme(  
        legend.title = element_blank(),
        axis.text = element_text(color = "black", size = 13),
        axis.text.x = element_text(size = 13),    # Hide axis text
        axis.title.y = element_text(size = 13),
        axis.ticks.y = element_line(color = "black"))+
        theme(axis.text = element_text( color = "black", size = 13), axis.title = element_text( size = 15 ))+
        theme(panel.background = element_rect(fill = "transparent", color = "black"))+ theme(plot.background = element_rect(color = NA))+
        theme(plot.title = element_text(hjust = 0.5, size = 15))

    # # calculate AUC
    bstype1 <- unique(roc_df$roc)[1]
    bstype2 <- unique(roc_df$roc)[2]
    auc1 <- round(roc_df$auc[roc_df$roc == bstype1][1],3)
    auc2 <- round(roc_df$auc[roc_df$roc == bstype2][2],3)

    p1 <- p1 + annotate("text", x = 0.8, y = 0.2, label = paste("AUC ", bstype1, ":", auc1), color = "blue", size = 4) +
        annotate("text", x = 0.8, y = 0.15, label = paste("AUC ", bstype2, ":", auc2), color = "red", size = 4)
    # save plot
    ggsave(path, plot=p1, width = 7, height = 7)
    return(p1)
}



# ----------------------------------------
#          References without IMD
# ----------------------------------------
tree_type = "ML"
bs_threshold = 80
bs_type = "ME"
ncol = 25
mode = "arithmetic_average"
references = c("ML")
# ROC df for ME+ML ref on bs_type ME and ME+IMD
roc_df_meimd <- prep_df_roc(tree_type = tree_type, ref = references, bs_threshold = bs_threshold, bs_type1 = "ME", bs_type2 = "ME+IMD", ncol = ncol, mode = mode)
p1 <- plot_ROC(roc_df_meimd, "plots/suppl/roc_me_meimd.png", title = "")



# ROC df for ME+ML ref on bs_type ML and ML+IMD
roc_df_mlmlimd <- prep_df_roc(tree_type = tree_type, ref = references, bs_threshold = bs_threshold, bs_type1 = "ML", bs_type2 = "ML+IMD", ncol = ncol, mode = mode)
p2 <- plot_ROC(roc_df_mlmlimd, "plots/suppl/roc_ml_mlmlimd.png", title = "")

# ROC df for ME+ML ref on bs_type ref ME+ML and ME+ML+IMD
roc_df_meiml <- prep_df_roc(tree_type = tree_type, ref = references, bs_threshold = bs_threshold, bs_type1 = "ME+ML", bs_type2 = "ME+ML+IMD", ncol = ncol, mode = mode)
p3 <- plot_ROC(roc_df_meiml, "plots/suppl/roc_meiml_meimlimd.png", title = "")


panel <- (p1+p2+p3)+ plot_annotation(tag_levels = 'A')+ plot_layout(ncol = 3)
# add title
ggsave(paste('plots/',"suppl",'/AUCs_refnoIMD.png', sep = ""), plot=panel, width = 18, height =6, dpi = 300)


# ----------------------------------------
#          References with IMD
# ----------------------------------------
tree_type = "ML"
bs_threshold = 80
bs_type = "ME"
ncol = 25
mode = "arithmetic_average"
references = c("IMD", "ME+IMD", "ML+IMD", "ME+ML+IMD")
# ROC df for ME+ML ref on bs_type ME and ME+IMD
roc_df_meimd <- prep_df_roc(tree_type = tree_type, ref = references, bs_threshold = bs_threshold, bs_type1 = "ME", bs_type2 = "ME+IMD", ncol = ncol, mode = mode)
p1 <- plot_ROC(roc_df_meimd, "plots/suppl/roc_me_meimd.png", title = "")



# ROC df for ME+ML ref on bs_type ML and ML+IMD
roc_df_mlmlimd <- prep_df_roc(tree_type = tree_type, ref = references, bs_threshold = bs_threshold, bs_type1 = "ML", bs_type2 = "ML+IMD", ncol = ncol, mode = mode)
p2 <- plot_ROC(roc_df_mlmlimd, "plots/suppl/roc_ml_mlmlimd.png", title = "")

# ROC df for ME+ML ref on bs_type ref ME+ML and ME+ML+IMD
roc_df_meiml <- prep_df_roc(tree_type = tree_type, ref = references, bs_threshold = bs_threshold, bs_type1 = "ME+ML", bs_type2 = "ME+ML+IMD", ncol = ncol, mode = mode)
p3 <- plot_ROC(roc_df_meiml, "plots/suppl/roc_meiml_meimlimd.png", title = "")


panel <- (p1+p2+p3)+ plot_annotation(tag_levels = 'A')+ plot_layout(ncol = 3)
# add title
ggsave(paste('plots/',"suppl",'/AUCs_refIMD.png', sep = ""), plot=panel, width = 18, height =6, dpi = 300)

# ----------------------------------------
#         All refs
# ----------------------------------------
tree_type = "ML"
bs_threshold = 80
ncol = 25
mode = "arithmetic_average"
references = c("ME", "ML", "ME+ML","IMD", "ME+IMD", "ML+IMD", "ME+ML+IMD")
# ROC df for ME+ML ref on bs_type ME and ME+IMD
roc_df_meimd <- prep_df_roc(tree_type = tree_type, ref = references, bs_threshold = bs_threshold, bs_type1 = "ME", bs_type2 = "ME+IMD", ncol = ncol, mode = mode)
p1 <- plot_ROC(roc_df_meimd, "plots/suppl/roc_me_meimd.png", title = "")

roc_df_mlmlimd

# ROC df for ME+ML ref on bs_type ML and ML+IMD
roc_df_mlmlimd <- prep_df_roc(tree_type = tree_type, ref = references, bs_threshold = bs_threshold, bs_type1 = "ML", bs_type2 = "ML+IMD", ncol = ncol, mode = mode)
p2 <- plot_ROC(roc_df_mlmlimd, "plots/suppl/roc_ml_mlmlimd.png", title = "")

# ROC df for ME+ML ref on bs_type ref ME+ML and ME+ML+IMD
roc_df_meiml <- prep_df_roc(tree_type = tree_type, ref = references, bs_threshold = bs_threshold, bs_type1 = "ME+ML", bs_type2 = "ME+ML+IMD", ncol = ncol, mode = mode)
p3 <- plot_ROC(roc_df_meiml, "plots/suppl/roc_meiml_meimlimd.png", title = "")


panel <- (p1+p2+p3)+ plot_annotation(tag_levels = 'A')+ plot_layout(ncol = 3)
# add title
ggsave(paste('plots/',"suppl",'/AUCs_allrefs.png', sep = ""), plot=panel, width = 18, height =6, dpi = 300)


roc_df_mlmlimd


# calculate auc per family 
tree_type = "ML"
bs_threshold = 80
ncol = 25
mode = "arithmetic_average"
ref = "ME+ML"
bs_type = "ME+IMD"
df_fams  <- data.frame()

for (fam in families){
    proven_positives = get_proven_positives(fam, tree_type, ref, bs_threshold)
    if(sum(proven_positives) == 0 | sum(proven_positives) == length(proven_positives)){
        next
    }
    predicted_values = get_splits_bs(fam, tree_type, bs_type, ncol,mode)                        
    auc_infos  = calculate_auc(proven_positives, predicted_values)
    auc_v = auc_infos$auc
    df_fams = rbind(df_fams, data.frame(fam = fam, auc = auc_v))
}

round(mean(df_fams$auc),3)
