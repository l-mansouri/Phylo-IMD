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
columns_filter = c("id", "dIMD", "dTM", "pML", "family", "percid", "seq1", "seq2")
df = input_patristic_df[,columns_filter]
df$pdist = 100-df$percid

# Normalize 
norm_imd = median(df$dIMD)
norm_tm = median(df$dTM)
norm_perc = median(df$pdist)
threshold = median(df$pML)
df$dIMD = df$dIMD*100/norm_imd
df$dTM= df$dTM*100/norm_tm
df$pdist = df$pdist*100/norm_perc
palette = c("#7e7e7e", "#A331D8", "#00ba39")
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
prep_tm_example
maxfill_imd <- max(hexbin(prep_imd_example$df[["pML"]], prep_imd_example$df[["dIMD"]], xbins = 50)@count)
maxfill_pdist <- max(hexbin(prep_perc_example$df[["pML"]], prep_perc_example$df[["pdist"]], xbins = 50)@count)
maxfill_tm <- max(hexbin(prep_tm_example$df[["pML"]], prep_tm_example$df[["dTM"]], xbins = 50)@count)
maxfillval = max(maxfill_pdist, maxfill_imd, maxfill_tm)

hexbin_example_imd <- plot_saturation_with_fitted_lines(prep_imd_example$df, prep_imd_example$line_before, prep_imd_example$line_after, "dIMD", "normalized IMD distance", threshold)+theme(legend.position = "right")
hexbin_example_tm <- plot_saturation_with_fitted_lines(prep_tm_example$df, prep_tm_example$line_before, prep_tm_example$line_after, "dTM", "normalized TM distance", threshold)+theme(legend.position = "right")
hexbin_example_perc <- plot_saturation_with_fitted_lines(prep_perc_example$df, prep_perc_example$line_before, prep_perc_example$line_after, "pdist", "normalized pdist", threshold)+theme(legend.position = "right")
hexbin_example_imd <- hexbin_example_imd+theme(legend.position = "right")+xlim(0,max(prep_tm_example$df$pML)+0.2)
hexbin_example_tm <-hexbin_example_tm+theme(legend.position = "none")+xlim(0,max(prep_tm_example$df$pML)+0.2)
hexbin_example_perc <- hexbin_example_perc+theme(legend.position = "none")+xlim(0,max(prep_tm_example$df$pML)+0.2)

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

  slope_threshold_imd <- round(coef(lm(as.formula(paste("dIMD", "~", "pML")), data = df_subset_threshold))[2],2)
  slope_threshold_tm <- round(coef(lm(as.formula(paste("dTM", "~", "pML")), data = df_subset_threshold))[2],2)
  slope_threshold_perc <- round(coef(lm(as.formula(paste("pdist", "~", "pML")), data = df_subset_threshold))[2],2)

  r2_imd <- round(cor(df_subset$dIMD, df_subset$pML)^2, 2)
  r2_tm <- round(cor(df_subset$dTM, df_subset$pML)^2, 2)
  r2_perc <- round(cor(df_subset$pdist, df_subset$pML)^2, 2)

  r2_imd_threshold <- round(cor(df_subset_threshold$dIMD, df_subset_threshold$pML)^2, 2)
  r2_tm_threshold <- round(cor(df_subset_threshold$dTM, df_subset_threshold$pML)^2, 2)
  r2_perc_threshold <- round(cor(df_subset_threshold$pdist, df_subset_threshold$pML)^2, 2)

  method = c("IMD", "TM", "pdist")
  slope = c(slope_imd, slope_tm, slope_perc)
  slope_threshold = c(slope_threshold_imd, slope_threshold_tm, slope_threshold_perc)
  r2 = c(r2_imd, r2_tm, r2_perc)
  r2_threshold = c(r2_imd_threshold, r2_tm_threshold, r2_perc_threshold)
  family = c(fam, fam, fam)
  slopes = data.frame(method, slope, r2,slope_threshold,r2_threshold, family, n_right_points, n_left_points)
  colnames(slopes) = c("method", "slope", "r2", "slope_threshold", "r2_threshold", "family", "n_right_points", "n_left_points")
  overall_slopes_per_family = rbind(overall_slopes_per_family, slopes)
}

overall_slopes_per_family$method = factor(overall_slopes_per_family$method,levels = c("pdist", "TM", "IMD"))
overall_slopes_per_family$ratio_slopes = (abs(overall_slopes_per_family$slope_threshold / overall_slopes_per_family$slope))

# -----------------------------------------------------------------------------
#  FILTERING
# -----------------------------------------------------------------------------
threshold_points = 20
overall_slopes_per_family = overall_slopes_per_family[overall_slopes_per_family$n_left_points > threshold_points,]
# right points
overall_slopes_per_family = overall_slopes_per_family[overall_slopes_per_family$n_right_points > threshold_points,]






# check how many families i have still
retained_families <- unique(overall_slopes_per_family$family)
length(retained_families)
#check how many sequences i have in total 
nseqs <- aggregate(id ~ family, data = df, FUN = length)
tot_seq <- sum(nseqs$id)
# total amount of sequences in retained families
sum(nseqs[nseqs$family %in% retained_families, "id"])*100/tot_seq


head(overall_slopes_per_family)
# save source data
write.table(overall_slopes_per_family, paste(source_data, "SATURATION_main_and_supp.csv", sep = ""), row.names = F, col.names = T)

head(overall_slopes_per_family)
# PLOT RATIOS
medians = aggregate(ratio_slopes ~ method, data = overall_slopes_per_family, FUN = median)
pratios <- ggplot(overall_slopes_per_family, aes(x = ratio_slopes, fill = method, color = method)) +
          geom_density(alpha = 0.5) +
          theme_minimal()+
          scale_fill_manual(values = palette)+
          scale_color_manual(values = palette)+
          geom_vline(data = medians, aes(xintercept = ratio_slopes, color = method), linetype = "dashed", linewidth = 1)+
          ggtitle("")+scale_x_log10()+xlab("Ratio slope before vs slope on all data (log10 scale)")+
          # remove legend title
          theme(legend.title = element_blank())+
          # all text blakc and size 10
          theme(axis.text = element_text(color = "black", size = 10))



c1 <-  p_sat1 + pratios + plot_annotation(tag_levels = list("A", "B"))
c1b <-  hexbin_example_perc + hexbin_example_tm + hexbin_example_imd + plot_annotation(tag_levels = list("C", "D", "E"))
combined_plot <- (c1b / c1) + plot_annotation(tag_levels = 'A')
combined_plot <- combined_plot+ plot_layout(heights = c(1.3, 2.2))
ggsave(paste(source_data, '../plots/review/', "SATURATION.png", sep = ''), plot = combined_plot, width = 11, height = 9, dpi = 300, units = 'in')


p_r2_before <- plot_distribution_across_families(overall_slopes_per_family, "r2_threshold")+xlab("R² before threshold")
p_r2_all <- plot_distribution_across_families(overall_slopes_per_family, "r2")+xlab("R²")
p_slope_before <- plot_distribution_across_families(overall_slopes_per_family, "slope_threshold")+xlab("Slope before threshold")
p_slope_all <- plot_distribution_across_families(overall_slopes_per_family, "slope")+xlab("Slope")


p1 <- (p_r2_before + p_r2_all)/(p_slope_before + p_slope_all)+ plot_annotation(tag_levels = 'A')
ggsave(paste(source_data, '../plots/review/', "SATURATION_SUPP.png", sep = ''), plot = p1, width = 10, height = 9, dpi = 300, units = 'in')


# make a table with the overall 



head(overall_slopes_per_family)
# remove family column, average by method all the other columns
overall_slopes_per_family_nofam = overall_slopes_per_family[, !(names(overall_slopes_per_family) %in% c("family", "n_right_points", "n_left_points"))]
table_data = aggregate(. ~ method, data = overall_slopes_per_family_nofam, FUN = mean)
# append here _mean to each column but method
colnames(table_data)[-1] = paste(colnames(table_data)[-1], "_mean", sep = "")
table_data_median = aggregate(. ~ method, data = overall_slopes_per_family_nofam, FUN = median)
colnames(table_data_median)[-1] = paste(colnames(table_data_median)[-1], "_median", sep = "")
table_data_sd = aggregate(. ~ method, data = overall_slopes_per_family_nofam, FUN = sd)
colnames(table_data_sd)[-1] = paste(colnames(table_data_sd)[-1], "_sd", sep = "")

table_data = merge(table_data, table_data_median, by = "method")
table_data = merge(table_data, table_data_sd, by = "method")

# order columns
colorder = c("method", "slope_threshold_mean", "slope_threshold_median", "slope_threshold_sd", "r2_threshold_mean", "r2_threshold_median", "r2_threshold_sd",  "slope_mean", "slope_median", "slope_sd", "r2_mean", "r2_median", "r2_sd")
table_data = table_data[,colorder]

# compute ratios on the median
table_data$ratio_slope = table_data$slope_threshold_median/table_data$slope_median
table_data$ratio_slope_median = table_data$slope_threshold_mean/table_data$slope_mean
# round everything but method by 2 decimals
table_data[,2:ncol(table_data)] = round(table_data[,2:ncol(table_data)], 2)


# save table
write.table(table_data, paste(source_data, "../tables/SATURATION_table.csv", sep = ""), row.names = F, col.names = T, sep = ",")
table_data
