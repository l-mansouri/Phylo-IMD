library(phangorn)
library(ggplot2)
library(ggExtra)
library(reshape2)
library(patchwork)
library(hexbin)
setwd('/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/')
source("04b_utils_saturation.R")

aligners=c('mTMalign')
trimming=c('untrimmed')
al=aligners[1]
tr=trimming[1]

source_data = '/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data/'
plots = paste(source_data, '../plots/', sep = '')
fl=read.table(paste(source_data, 'list_of_families_with_all_rep_in_3d', sep = '/'))[,1]

input_patristic_df = read.table(paste(source_data,al,'_',tr,'_input_patristic_distances_test.txt', sep=''), header = TRUE)
perc_sim = read.table(paste(source_data, 'percids.txt', sep = ''), header = TRUE)

# -----------------------------------------------------------------------------
#  PREPARE DATAFRAMES TO BE MERGED
# -----------------------------------------------------------------------------
# PREPARE PERC IDENTITY DATA FRAME
# create a new column with seqquence 1 and 2 ordered in alphaebtical order
perc_sim$id = pmin(perc_sim$seq1, perc_sim$seq2)
perc_sim$id = paste(perc_sim$id, pmax(perc_sim$seq1, perc_sim$seq2), sep = "_")
perc_sim$id = paste(perc_sim$id, perc_sim$family, sep = "_")
perc_sim$id = gsub(".pdb", ".p", perc_sim$id)
perc_sim = perc_sim[, !(names(perc_sim) %in% c("family"))]

# now do the same for the input_patristic_df
input_patristic_df$seq1 = sapply(strsplit(input_patristic_df$seqs, "_"), "[", 1)
input_patristic_df$seq2 = sapply(strsplit(input_patristic_df$seqs, "_"), "[", 2)
input_patristic_df$id = pmin(input_patristic_df$seq1, input_patristic_df$seq2)
input_patristic_df$id = paste(input_patristic_df$id, pmax(input_patristic_df$seq1, input_patristic_df$seq2), sep = "_")
input_patristic_df$id = paste(input_patristic_df$id, input_patristic_df$family, sep = "_")
input_patristic_df = input_patristic_df[, !(names(input_patristic_df) %in% c("seq1", "seq2"))]

# now merge the two data frames
input_patristic_df = merge(input_patristic_df, perc_sim, by = "id", all.x = TRUE)

# Calculate the slope before and after the 0.25 threshold om the ML patristic distance
columns_filter = c("id", "dIMD", "dTM", "pML", "family", "percid", "dME", "seq1", "seq2")
df = input_patristic_df[,columns_filter]
df$pdist = 100-df$percid

# Normalize 
norm_imd = median(df$dIMD)
norm_tm = median(df$dTM)
norm_perc = median(df$pdist)
threshold = median(df$pML)
norm_me = median(df$dME)
df$dIMD = df$dIMD*100/norm_imd
df$dTM= df$dTM*100/norm_tm
df$pdist = df$pdist*100/norm_perc
df$dME = df$dME*100/norm_me
palette = c("#7e7e7e", "#f8766d", "#A331D8", "#00ba39")
color_line = "#4ca0ff"
pseudocount = 0
threshold_filtering = 0.3
line_size = 1.3




# -----------------------------------------------------------------------------
#       PLOT example
# -----------------------------------------------------------------------------
example_fam = "PF13378"
df_example = df[df$family == example_fam,]
prep_imd_example <- get_before_after(df_example, "dIMD", threshold)
prep_tm_example <- get_before_after(df_example, "dTM", threshold)
prep_perc_example <- get_before_after(df_example, "pdist", threshold)
prep_me_example <- get_before_after(df_example, "dME", threshold)
prep_tm_example
maxfill_imd <- max(hexbin(prep_imd_example$df[["pML"]], prep_imd_example$df[["dIMD"]], xbins = 50)@count)
maxfill_pdist <- max(hexbin(prep_perc_example$df[["pML"]], prep_perc_example$df[["pdist"]], xbins = 50)@count)
maxfill_tm <- max(hexbin(prep_tm_example$df[["pML"]], prep_tm_example$df[["dTM"]], xbins = 50)@count)
maxfill_me <- max(hexbin(prep_me_example$df[["pML"]], prep_me_example$df[["dME"]], xbins = 50)@count)
maxfillval = max(maxfill_pdist, maxfill_imd, maxfill_tm, maxfill_me)

hexbin_example_imd <- plot_saturation_with_fitted_lines(prep_imd_example$df, prep_imd_example$line_before, prep_imd_example$line_after, "dIMD", "normalized IMD distance", threshold)+theme(legend.position = "right")
hexbin_example_tm <- plot_saturation_with_fitted_lines(prep_tm_example$df, prep_tm_example$line_before, prep_tm_example$line_after, "dTM", "normalized TM distance", threshold)+theme(legend.position = "right")
hexbin_example_perc <- plot_saturation_with_fitted_lines(prep_perc_example$df, prep_perc_example$line_before, prep_perc_example$line_after, "pdist", "pdist", threshold)+theme(legend.position = "right")
hexbin_example_me <- plot_saturation_with_fitted_lines(prep_me_example$df, prep_me_example$line_before, prep_me_example$line_after, "dME", "normalized ME distance", threshold)+theme(legend.position = "right")
hexbin_example_imd <- hexbin_example_imd+theme(legend.position = "right")+xlim(0,max(prep_tm_example$df$pML)+0.2)
hexbin_example_tm <-hexbin_example_tm+theme(legend.position = "none")+xlim(0,max(prep_tm_example$df$pML)+0.2)
hexbin_example_perc <- hexbin_example_perc+theme(legend.position = "none")+xlim(0,max(prep_tm_example$df$pML)+0.2)
hexbin_example_me <- hexbin_example_me+theme(legend.position = "none")+xlim(0,max(prep_tm_example$df$pML)+0.2)

# source data
write.table(df_example, paste(source_data, "SATURATION_example_main.csv", sep = ""), row.names = F, col.names = T)




# -----------------------------------------------------------------------------
#       PLOT SLOPES
# -----------------------------------------------------------------------------


overall_slopes_per_family = data.frame()
for (fam in fl){
  df_subset = df[df$family == fam,]
  # df up to threshold
  df_subset_threshold = df_subset[df_subset$pML < threshold,]

  n_right_points = nrow(df_subset) - nrow(df_subset_threshold)
  n_left_points = nrow(df_subset_threshold)



  # get slope for imd vs pML
  slope_imd <- round(coef(lm(as.formula(paste("dIMD", "~", "pML")), data = df_subset))[2],2)
  slope_tm <- round(coef(lm(as.formula(paste("dTM", "~", "pML")), data = df_subset))[2],2)
  slope_perc <- round(coef(lm(as.formula(paste("pdist", "~", "pML")), data = df_subset))[2],2)
  slope_me <- round(coef(lm(as.formula(paste("dME", "~", "pML")), data = df_subset))[2],2)

  slope_threshold_imd <- round(coef(lm(as.formula(paste("dIMD", "~", "pML")), data = df_subset_threshold))[2],2)
  slope_threshold_tm <- round(coef(lm(as.formula(paste("dTM", "~", "pML")), data = df_subset_threshold))[2],2)
  slope_threshold_perc <- round(coef(lm(as.formula(paste("pdist", "~", "pML")), data = df_subset_threshold))[2],2)
  slope_threshold_me <- round(coef(lm(as.formula(paste("dME", "~", "pML")), data = df_subset_threshold))[2],2)

  r2_imd <- round(cor(df_subset$dIMD, df_subset$pML)^2, 2)
  r2_tm <- round(cor(df_subset$dTM, df_subset$pML)^2, 2)
  r2_perc <- round(cor(df_subset$pdist, df_subset$pML)^2, 2)
  r2_me <- round(cor(df_subset$dME, df_subset$pML)^2, 2)

  r2_imd_threshold <- round(cor(df_subset_threshold$dIMD, df_subset_threshold$pML)^2, 2)
  r2_tm_threshold <- round(cor(df_subset_threshold$dTM, df_subset_threshold$pML)^2, 2)
  r2_perc_threshold <- round(cor(df_subset_threshold$pdist, df_subset_threshold$pML)^2, 2)
  r2_me_threshold <- round(cor(df_subset_threshold$dME, df_subset_threshold$pML)^2, 2)

  method = c("IMD", "TM", "pdist", "ME")
  slope = c(slope_imd, slope_tm, slope_perc, slope_me)
  slope_threshold = c(slope_threshold_imd, slope_threshold_tm, slope_threshold_perc, slope_threshold_me)
  r2 = c(r2_imd, r2_tm, r2_perc, r2_me)
  r2_threshold = c(r2_imd_threshold, r2_tm_threshold, r2_perc_threshold, r2_me_threshold)
  n_points_interval = c(n_points_interval_imd, n_points_interval_tm, n_points_interval_perc, n_points_interval_me)
  family = c(fam, fam, fam, fam)
  slopes = data.frame(method, slope, r2,slope_threshold,r2_threshold, family, n_right_points, n_left_points, n_points_interval)
  colnames(slopes) = c("method", "slope", "r2", "slope_threshold", "r2_threshold", "family", "n_right_points", "n_left_points")
  overall_slopes_per_family = rbind(overall_slopes_per_family, slopes)
}

overall_slopes_per_family$method = factor(overall_slopes_per_family$method,levels = c("pdist", "ME", "TM", "IMD"))
overall_slopes_per_family$ratio_slopes = (abs(overall_slopes_per_family$slope_threshold / overall_slopes_per_family$slope))

# -----------------------------------------------------------------------------
#  FILTERING
# -----------------------------------------------------------------------------
threshold_points = 20
overall_slopes_per_family = overall_slopes_per_family[overall_slopes_per_family$n_left_points > threshold_points,]
overall_slopes_per_family = overall_slopes_per_family[overall_slopes_per_family$n_right_points > threshold_points,]


head(overall_slopes_per_family)

# check how many families i have still
retained_families <- unique(overall_slopes_per_family$family)
length(retained_families)
#check how many sequences i have in total 
nseqs <- aggregate(id ~ family, data = df, FUN = length)
tot_seq <- sum(nseqs$id)
# total amount of sequences in retained families
sum(nseqs[nseqs$family %in% retained_families, "id"])*100/tot_seq

# save source data
write.table(overall_slopes_per_family, paste(source_data, "SATURATION_main_and_supp.csv", sep = ""), row.names = F, col.names = T)

# rename in method pdist to ME
# factor methods
overall_slopes_per_family$method = factor(overall_slopes_per_family$method,levels = c("pdist", "ME", "TM", "IMD"))
# PLOT RATIOS
medians = aggregate(ratio_slopes ~ method, data = overall_slopes_per_family, FUN = median)
pratios <- ggplot(overall_slopes_per_family, aes(x = ratio_slopes, fill = method, color = method)) +
          geom_density(alpha = 0.5) +
          theme_minimal()+
          scale_fill_manual(values = palette)+
          scale_color_manual(values = palette)+
          geom_vline(data = medians, aes(xintercept = ratio_slopes, color = method), linetype = "dashed", linewidth = 1)+
          ggtitle("")+scale_x_log10()+xlab("Ratio of linear regression slopes (homologues/all) \n (log10 scale)")+
          # remove legend title
          theme(legend.title = element_blank())+
          # all text blakc and size 10
          theme(axis.text = element_text(color = "black", size = 10))
# rename ME to pdist

# get a new df with each value in column method a new column

df_for_test = dcast(overall_slopes_per_family, family ~ method, value.var = "ratio_slopes")

# paired wilcoxon test on the ratios btw TM and pdist
wilcox.test(df_for_test$TM, df_for_test$pdist, paired = TRUE)

# paired wilcoxon test on the ratios btw IMD and pdist
wilcox.test(df_for_test$IMD, df_for_test$pdist, paired = TRUE)



c1 <-  pratios + p_sat1 + plot_layout(widths = c(1.8, 1))
hexbins_sequence <-  hexbin_example_perc + hexbin_example_me
hexbins_structures <-  hexbin_example_tm + hexbin_example_imd
combined_plot <- (hexbins_sequence / hexbins_structures / c1) + plot_annotation(tag_levels = 'A')
combined_plot <- combined_plot+ plot_layout(heights = c(1.3,1.3, 1.5))
ggsave(paste(source_data, '../plots/review/', "SATURATION.png", sep = ''), plot = combined_plot, width = 11, height = 15, dpi = 300, units = 'in')


p_r2_before <- plot_distribution_across_families(overall_slopes_per_family, "r2_threshold")+xlab("R² (close homologues)")
p_r2_all <- plot_distribution_across_families(overall_slopes_per_family, "r2")+xlab("R²")
p_slope_before <- plot_distribution_across_families(overall_slopes_per_family, "slope_threshold")+xlab("Slope (close homologues)")
p_slope_all <- plot_distribution_across_families(overall_slopes_per_family, "slope")+xlab("Slope")


p1 <- (p_slope_before + p_slope_all)/(p_r2_before + p_r2_all)+ plot_annotation(tag_levels = 'A')
ggsave(paste(source_data, '../plots/review/', "SATURATION_SUPP.png", sep = ''), plot = p1, width = 10, height = 9, dpi = 300, units = 'in')
# save source data
write.table(overall_slopes_per_family, paste(source_data, "SATURATION_supp.csv", sep = ""), row.names = F, col.names = T, sep = ",")



# -----------------------------------------------------------------------------
#  CREATE TABLE
# -----------------------------------------------------------------------------

# remove family column, average by method all the other columns
head(overall_slopes_per_family)
# how many ratio slopes between 0.8 and 1.2, separated by method
total_families_considered = length(unique(overall_slopes_per_family$family))
overall_slopes_per_family$within_range = ifelse(overall_slopes_per_family$ratio_slopes > 0.8 & overall_slopes_per_family$ratio_slopes < 1.2, 1, 0)
overall_slopes_per_family_for_sum  = overall_slopes_per_family[,c("method", "within_range")]
table_within_range = aggregate(. ~ method, data = overall_slopes_per_family_for_sum, FUN = sum)
table_within_range$within_range = table_within_range$within_range*100/total_families_considered
table_within_range


overall_slopes_per_family_nofam = overall_slopes_per_family[, !(names(overall_slopes_per_family) %in% c("family", "n_right_points", "n_left_points"))]
table_data_median = aggregate(. ~ method, data = overall_slopes_per_family_nofam, FUN = median)
colnames(table_data_median)[-1] = paste(colnames(table_data_median)[-1], "_median", sep = "")
table_data_sd = aggregate(. ~ method, data = overall_slopes_per_family_nofam, FUN = sd)
colnames(table_data_sd)[-1] = paste(colnames(table_data_sd)[-1], "_sd", sep = "")

table_data = merge(table_within_range, table_data_median, by = "method")
table_data = merge(table_data, table_data_sd, by = "method")
head(table_data)

# order columns
colorder = c("method", "slope_threshold_median", "slope_threshold_sd", "r2_threshold_median", "r2_threshold_sd",   "slope_median", "slope_sd", "r2_median", "r2_sd", "ratio_slopes_median", "ratio_slopes_sd", "within_range")
table_data = table_data[,colorder]

# round everything but method by 2 decimals
table_data[,2:ncol(table_data)] = round(table_data[,2:ncol(table_data)], 2)
# order rows table data by method (pdist, ME, TM, IMD)
table_data = table_data[order(table_data$method),]
table_data

# save table
write.table(table_data, paste(source_data, "../tables/SATURATION_table.csv", sep = ""), row.names = F, col.names = T, sep = ",")



# --------------------------------------------------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#  GET EXAMPLES
# -----------------------------------------------------------------------------
plot_scatters_distances <- function(df, low_5, plot_name){
  # plot the low families
  plot_low_5_list = list()
  i = 1 
  for (fam in low_5){

    
    # compute R
    input_patristic_df_subset = df[df$family == fam,]


    # add perc id
    r = cor(input_patristic_df_subset$pdist, input_patristic_df_subset$pML, method = "pearson")
    r2 = r^2
    plot_sat = saturation_plot_smooth(input_patristic_df_subset, "pML", "pdist", "ML patristic distance", "pdist")+ggtitle(paste("R² = ", round(r2, 2)))
    plot_sat = plot_sat +ggtitle(paste(fam, "\n R² = ", round(r2, 2), "\n pdist = ", round(mean(input_patristic_df_subset$pdist), 2)))
    # put title in the milddle
    plot_sat = plot_sat + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
    plot_low_5_list[[i]] = plot_sat

    # ME
    r = cor(input_patristic_df_subset$dME, input_patristic_df_subset$pML, method = "pearson")
    r2 = r^2
    plot_sat = saturation_plot_smooth(input_patristic_df_subset, "pML", "dME", "ML patristic distance", "normalized ME distance")+ggtitle(paste("R² = ", round(r2, 2)))
    plot_sat = plot_sat +ggtitle(paste(fam, "\n R² = ", round(r2, 2), "\n pdist = ", round(mean(input_patristic_df_subset$pdist), 2)))
    # put title in the milddle
    plot_sat = plot_sat + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
    plot_low_5_list[[i+1]] = plot_sat


    r = cor(input_patristic_df_subset$dTM, input_patristic_df_subset$pML, method = "pearson")
    r2 = r^2
    plot_sat = saturation_plot_smooth(input_patristic_df_subset, "pML", "dTM", "ML patristic distance", "normalized TM distance")+ggtitle(paste("R² = ", round(r2, 2)))
    plot_sat = plot_sat +ggtitle(paste(fam, "\n R² = ", round(r2, 2), "\n pdist = ", round(mean(input_patristic_df_subset$pdist), 2)))
    # put title in the milddle
    plot_sat = plot_sat + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
    plot_low_5_list[[i+2]] = plot_sat


    r = cor(input_patristic_df_subset$dIMD, input_patristic_df_subset$pML, method = "pearson")
    r2 = r^2
    plot_sat = saturation_plot_smooth(input_patristic_df_subset, "pML", "dIMD", "ML patristic distance", "normalized IMD distance")
    plot_sat = plot_sat +ggtitle(paste(fam, "\n R² = ", round(r2, 2), "\n pdist = ", round(mean(input_patristic_df_subset$pdist), 2)))
    plot_sat = plot_sat + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
    plot_low_5_list[[i+3]] = plot_sat
    i = i + 4


  }

  # put the low 5 in a panel with 2 col using patchwork
  # Combine plots into a grid with one common legend
  combined_plot <- (plot_low_5_list[[1]]+
                    plot_low_5_list[[2]] + 
                    plot_low_5_list[[3]] +
                    plot_low_5_list[[4]] +
                    plot_low_5_list[[5]] +
                    plot_low_5_list[[6]] +
                    plot_low_5_list[[7]] +
                    plot_low_5_list[[8]] +
                    plot_low_5_list[[9]] +
                    plot_low_5_list[[10]]+
                    plot_low_5_list[[11]]+
                    plot_low_5_list[[12]]+
                    plot_low_5_list[[13]]+
                    plot_low_5_list[[14]]+
                    plot_low_5_list[[15]]+
                    plot_low_5_list[[16]]+
                    plot_low_5_list[[17]]+
                    plot_low_5_list[[18]]+
                    plot_low_5_list[[19]]+
                    plot_low_5_list[[20]]
                    ) + 
    plot_layout(ncol = 4, guides = "collect") & 
    theme(legend.position = 'right')&
    plot_annotation(tag_levels = 'A')

  # save the final image
  ggsave(paste(source_data, '../plots/review/', plot_name, sep = ''), plot = combined_plot, width = 12, height = 18, dpi = 300, units = 'in')

}

df = input_patristic_df[,columns_filter]
df$pdist = 100-df$percid

# Normalize 
norm_imd = median(df$dIMD)
norm_tm = median(df$dTM)
norm_perc = median(df$pdist)
norm_me = median(df$dME)
threshold = median(df$pML)
df$dIMD = df$dIMD*100/norm_imd
df$dTM= df$dTM*100/norm_tm
df$pdist = df$pdist*100/norm_perc
df$dME = df$dME*100/norm_me


correlations = data.frame()
for (fam in fl){

  df_subset = df[df$family == fam,]
  cor_dIMD_pML = cor(df_subset$dIMD, df_subset$pML, method = "pearson")
  r2_imd_ml = cor_dIMD_pML^2

  cor_dTM_pML = cor(df_subset$dTM, df_subset$pML, method = "pearson")
  r2_tm_ml = cor_dTM_pML^2

  cor_percid_pML = cor(df_subset$pdist, df_subset$pML, method = "pearson")
  r2_percid_ml = cor_percid_pML^2

  correlations = rbind(correlations, data.frame(family = fam, cor_dIMD_pML = cor_dIMD_pML, r2_imd_ml = r2_imd_ml, cor_dTM_pML = cor_dTM_pML, r2_tm_ml = r2_tm_ml, cor_percid_pML = cor_percid_pML, r2_percid_ml = r2_percid_ml))
}



# Get bottom at 9th decile
perc_quantile_threshold_10 = quantile(correlations$cor_dIMD_pML, c(0.1)) 
removed_quantile_10 = correlations[correlations$cor_dIMD_pML > perc_quantile_threshold_10,]
removed_quantile_10 = removed_quantile_10[order(removed_quantile_10$cor_dIMD_pML),]
low_5 = removed_quantile_10[1:5,"family"]
low_5

# Get middle at 5th decile
perc_quantile_threshold_50 = quantile(correlations$cor_dIMD_pML, c(0.5, 0.4))
# get the values in the middle
quantile_50_elements = correlations[correlations$cor_dIMD_pML > perc_quantile_threshold_50[2] & correlations$cor_dIMD_pML < perc_quantile_threshold_50[1],]
# now get the middle ones, sort and get the middle
# order 
quantile_50_elements = quantile_50_elements[order(quantile_50_elements$cor_dIMD_pML),]
middle_5th_idx = round(nrow(quantile_50_elements)/2)
lower_middle_5 = middle_5th_idx - 2
upper_middle_5 = middle_5th_idx + 2
mid_5 = quantile_50_elements[lower_middle_5:upper_middle_5,"family"]
mid_5

# get 5 random elements near the median 
# extract the median
med = median(correlations$cor_dIMD_pML)
# get the 5 closest values to the median
correlations$diff = abs(correlations$cor_dIMD_pML - med)
# sort by diff
correlations = correlations[order(correlations$diff),]
mid_5 = correlations[1:5,]$family

# now extract the 5 families with the most number of sequences
# get the number of sequences per family
seqs_per_fam = data.frame()
for (fam in fl){
  seqs_per_fam = rbind(seqs_per_fam, data.frame(family = fam, n_seqs = nrow(input_patristic_df[input_patristic_df$family == fam,])))
}

# get the 5 families with the most sequences
top_5_powered = seqs_per_fam[order(seqs_per_fam$n_seqs, decreasing = TRUE),][1:5,"family"]
top_5_powered

# get the ones with the highest correlation
top_5 = correlations[order(correlations$cor_dIMD_pML, decreasing = TRUE),][1:5,"family"]



#--------------------------------------------------------------------------
# PLOT THE SCATTERS
#--------------------------------------------------------------------------
plot_scatters_distances(df, low_5, 'SATURATION_bottom_9th_decile.png')
plot_scatters_distances(df, mid_5, 'SATURATION_middle_5th_decile.png')
plot_scatters_distances(df, top_5_powered, 'SATURATION_top_5_powered.png')
plot_scatters_distances(df, top_5, 'SATURATION_top_5.png')

#--------------------------------------------------------------------------
# SAVE SOURCE DATA
#--------------------------------------------------------------------------
input_patristic_df_for_sd = df[,c("family", "seq1","seq2", "dIMD", "dTM", "dME" ,"pdist", "pML")]

head(df)

# save the input patristic distances for low_5 families in a file
write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% low_5,], paste(source_data, 'saturation_low_5.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")
# save the input patristic distances for mid_5 families in a file
write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% mid_5,], paste(source_data, 'saturation_mid_5.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")
# save the input patristic distances for top_5 families in a file
write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% top_5,], paste(source_data, 'saturation_top_5.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")
# save the input patristic distances for top_5 families in a file
write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% top_5_powered,], paste(source_data, 'saturation_top_5_powered.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")
