library(phangorn)
library(ggplot2)
library(ggExtra)
library(reshape2)
library(patchwork)
library(hexbin)


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

# -----------------------------------------------------------------------------
#  ACCESSORY PLOTTING FUNCTIONS
# -----------------------------------------------------------------------------
saturation_plot <- function(input_patristic_df, distance, patristic, distance_label = "", patristic_label = "", maxfillval = NA){
  plot_saturation =ggplot(input_patristic_df, aes_string(x=distance, y=patristic))+
          geom_point(alpha=0)+
          geom_hex(bins=50)+
          xlab(paste(distance_label, '', sep = " "))+
          ylab(paste(patristic_label, '', sep = " "))+
          theme_light()+
          scale_fill_gradientn(colors = c("#d4d4d4", rgb(0, 0.05, 0.3)), trans = "log", limits = c(NA, maxfillval),labels = function(x) format(round(x, 0), scientific = FALSE)) +
          theme(axis.text = element_text(color = "black", size = 10),
                axis.text.x = element_text(size = 10),
                axis.title.y = element_text(size = 15),
                axis.title.x = element_text(size = 15),
                axis.ticks.y = element_line(color = "black"))+
          expand_limits(y=0)+
          expand_limits(x=0)

  #ggsave(paste(plots,al,'_',tr,'_',x_label,'_', y_label,'_saturation.png', sep=''), plot=plot_saturation)
  return(plot_saturation)
}

saturation_plot_smooth <- function(input_patristic_df, distance, patristic, distance_label = "", patristic_label = ""){
  plot_saturation =ggplot(input_patristic_df, aes_string(x=distance, y=patristic))+
          geom_point(alpha=0.3)+
          geom_smooth(method = "lm")+
          xlab(paste(distance_label, '', sep = " "))+
          ylab(paste(patristic_label, '', sep = " "))+
          theme_light()+
          theme(axis.text = element_text(color = "black", size = 10),
                axis.text.x = element_text(size = 10),
                axis.title.y = element_text(size = 15),
                axis.title.x = element_text(size = 15),
                axis.ticks.y = element_line(color = "black"))+
          expand_limits(y=0)+
          expand_limits(x=0)
  #ggsave(paste(plots,al,'_',tr,'_',distance_label,'_', patristic_label,'_saturation_smooth.png', sep=''), plot=plot_saturation)
  return(plot_saturation)
}


get_before_after <- function(df, y, threshold){
  # IMD vs ML patristic distance
  df_imd_ml = df[,c(y, "pML", "family", "seq1", "seq2")]
  df_imd_ml_before = df_imd_ml[df_imd_ml$pML < threshold,]
  df_imd_ml_after = df_imd_ml[df_imd_ml$pML >= threshold,]
  # Check if there are enough points to fit a line
  if(nrow(df_imd_ml_before) < 2 | nrow(df_imd_ml_after) < 2){
    return( list(df = df_imd_ml, line_before = NA, line_after = NA, slope_before = NA, slope_after = NA, r2_before = NA, r2_after = NA))
  }
  line_imd_ml_before = lm(as.formula(paste(y, "~", "pML")), data = df_imd_ml_before)
  slope_imd_ml_before <- coef(line_imd_ml_before)[2]
  line_imd_ml_after = lm(as.formula(paste(y, "~", "pML")), data = df_imd_ml_after)
  slope_imd_ml_after <- coef(line_imd_ml_after)[2]
  # round
  slope_imd_ml_before = round(slope_imd_ml_before, 2)
  slope_imd_ml_after = round(slope_imd_ml_after, 2)
  # calculate R²
  r_imd_ml_before = cor(df_imd_ml_before[[y]], df_imd_ml_before$pML, method = "pearson")
  r2_imd_ml_before = round(r_imd_ml_before^2,2)
  r_imd_ml_after = cor(df_imd_ml_after[[y]], df_imd_ml_after$pML, method = "pearson")
  r2_imd_ml_after = round(r_imd_ml_after^2,2)

  return( list(df = df_imd_ml, line_before = line_imd_ml_before, line_after = line_imd_ml_after, slope_before = slope_imd_ml_before, slope_after = slope_imd_ml_after, r2_before = r2_imd_ml_before, r2_after = r2_imd_ml_after))

}



# Calculate the slope before and after the 0.25 threshold om the ML patristic distance
columns_filter = c("id", "dIMD", "dTM", "pML", "family", "percid", "seq1", "seq2")
df = input_patristic_df[,columns_filter]
df$pdist = 100-df$percid

# Normalize 
norm_imd = median(df$dIMD)
norm_tm = median(df$dTM)
norm_perc = median(df$pdist)
threshold = median(df$pML)
df$dIMD = df$dIMD/norm_imd
df$dTM= df$dTM/norm_tm
df$pdist = df$pdist/norm_perc
palette = c("#7e7e7e","#ffaf1c", "#2191FB" )
palette = c("#7e7e7e", "#A331D8", "#00ba39")
color_line = "#41e9ff"
pseudocount = 0.001

slopes_per_family = data.frame()
for (fam in fl){
  df_subset = df[df$family == fam,]
  imd_subset_prep = get_before_after(df_subset, "dIMD", threshold)

  tm_subset_prep = get_before_after(df_subset, "dTM", threshold)

  perc_subset_prep = get_before_after(df_subset, "pdist", threshold)

  slopes = data.frame(method = c("IMD ", "TM ", "pdist"),
                   before = c(imd_subset_prep$slope_before, tm_subset_prep$slope_before, perc_subset_prep$slope_before),
                   after = c(imd_subset_prep$slope_after, tm_subset_prep$slope_after, perc_subset_prep$slope_after),
                    r2_before = c(imd_subset_prep$r2_before, tm_subset_prep$r2_before, perc_subset_prep$r2_before),
                    r2_after = c(imd_subset_prep$r2_after, tm_subset_prep$r2_after, perc_subset_prep$r2_after))

  slopes$ratio = round(slopes$before+pseudocount/slopes$after+pseudocount, 2)
  slopes$diff = slopes$before - slopes$after
  slopes$family = fam
  slopes_per_family = rbind(slopes_per_family, slopes)
}

# Check how many families i have removed from the analysis
families_excluded <- (unique(slopes_per_family$family[is.na(slopes_per_family$before)]))
families_excluded_after <- (unique(slopes_per_family$family[is.na(slopes_per_family$after)]))
slopes_per_family = slopes_per_family[!is.na(slopes_per_family$before),]
slopes_per_family = slopes_per_family[!is.na(slopes_per_family$after),]

plot_distribution_across_families <- function(slopes_per_family, y){
  correlations_melted = melt(slopes_per_family[,c(y, "method")])
  # plot the distribution of R² values before
  correlations_melted$method = factor(correlations_melted$method,levels = c("pdist", "TM ", "IMD "))
  slopes_per_family$y = slopes_per_family[[y]]
  medians = aggregate(y ~ method, data = slopes_per_family, FUN = median)
  p <- ggplot(correlations_melted, aes(x = value, fill = method, color = method)) +
              geom_density(alpha = 0.3) +
              theme_light() +
              xlab(y) + 
              ylab("Density") +
              theme(legend.title=element_blank())+ theme(legend.position = "right")+
              geom_vline(data = medians, aes(xintercept = y, color = method), linetype = "dashed", linewidth = 1)+
              scale_fill_manual(values = palette)+
              scale_color_manual(values = palette)
  return(p)
}



plot_boxplot_across_families <- function(slopes_per_family, y){
  correlations_melted = melt(slopes_per_family[,c(y, "method")])
  # plot the distribution of R² values before
  correlations_melted$method = factor(correlations_melted$method, levels = c("pdist", "TM ", "IMD "))
  slopes_per_family$y = slopes_per_family[[y]]
  p <- ggplot(correlations_melted, aes(x = method, y = value, fill = method, color = method)) +
              geom_boxplot(alpha = 0.5) +
              theme_light() +
              xlab("") + 
              ylab(paste(y)) +
              theme(legend.title=element_blank())+ theme(legend.position = "right")+
              scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
              scale_fill_manual(values = palette)+
              scale_color_manual(values = palette)+
              scale_y_continuous(labels = function(x) format(x, scientific = FALSE))
  return(p)
}


# -----------------------------------------------------------------------------
#       PERC ID
# -----------------------------------------------------------------------------


line_size = 1.3
plot_saturation_with_fitted_lines <- function(df_imd_ml,line_imd_ml_before, line_imd_ml_after, y, ylabel,  threshold, maxfillval = NA){
  threshold_x <- threshold
  IMD_intercept <- coef(line_imd_ml_before)[1]
  IMD_slope <- coef(line_imd_ml_before)[2]
  IMD_threshold_y <- IMD_intercept + IMD_slope * threshold_x
  IMD_max_y = max(df_imd_ml[[y]])
  IMD_x_for_y_100 = (IMD_max_y - IMD_intercept)/IMD_slope
  

  IMD_intercept_after <- coef(line_imd_ml_after)[1]
  IMD_slope_after <- coef(line_imd_ml_after)[2]
  IMD_threshold_y_after <- IMD_intercept_after + IMD_slope_after * threshold_x

  pml_dimd <- saturation_plot(df_imd_ml, "pML", y, "ML patristic distance", ylabel, maxfillval) +
    # Line before the threshold (solid)
    geom_segment(aes(x = min(df_imd_ml$pML), y = IMD_intercept + IMD_slope * min(df_imd_ml$pML), 
                    xend = threshold_x, yend = IMD_threshold_y), color = color_line, size =line_size) +
    # Line after the threshold (dashed)
    geom_segment(aes(x = threshold_x, y = IMD_threshold_y, 
                    xend = IMD_x_for_y_100, yend = IMD_max_y), 
                color = color_line, size =line_size, linetype = "dashed")+
    geom_segment(aes(x = min(df_imd_ml$pML), y = IMD_intercept_after + IMD_slope_after * min(df_imd_ml$pML), 
                    xend = threshold_x, yend = IMD_threshold_y_after), color = color_line, size =line_size, linetype = "dashed") +
    # Line after the threshold (dashed)
    geom_segment(aes(x = threshold_x, y = IMD_threshold_y_after, 
                    xend = max(df_imd_ml$pML), yend = IMD_intercept_after + IMD_slope_after * max(df_imd_ml$pML)), 
                color = color_line, size =line_size)

  r = cor(df_imd_ml[[y]], df_imd_ml$pML, method = "pearson")
  r2 = r^2
  pml_dimd = pml_dimd + ggtitle(paste("R² = ", round(r2, 2)))
  pml_dimd = pml_dimd + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
  pml_dimd = pml_dimd + geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 0.7)
  return(pml_dimd)

}


# -----------------------------------------------------------------------------
#       PLOT example
# -----------------------------------------------------------------------------
example_fam = "PF13378"
df_example = df[df$family == example_fam,]
prep_imd_example <- get_before_after(df_example, "dIMD", threshold)
prep_tm_example <- get_before_after(df_example, "dTM", threshold)
prep_perc_example <- get_before_after(df_example, "pdist", threshold)

maxfill_imd <- max(hexbin(prep_imd_example$df[["pML"]], prep_imd_example$df[["dIMD_norm"]], xbins = 50)@count)
maxfill_pdist <- max(hexbin(prep_perc_example$df[["pML"]], prep_perc_example$df[["pdist_norm"]], xbins = 50)@count)
maxfill_tm <- max(hexbin(prep_tm_example$df[["pML"]], prep_tm_example$df[["dTM_norm"]], xbins = 50)@count)
maxfillval = max(maxfill_pdist, maxfill_imd, maxfill_tm)
maxfillval
hexbin_example_imd <- plot_saturation_with_fitted_lines(prep_imd_example$df, prep_imd_example$line_before, prep_imd_example$line_after, "dIMD", "normalized IMD distance", threshold)+theme(legend.position = "right")
hexbin_example_tm <- plot_saturation_with_fitted_lines(prep_tm_example$df, prep_tm_example$line_before, prep_tm_example$line_after, "dTM", "normalized TM distance", threshold)+theme(legend.position = "right")
hexbin_example_perc <- plot_saturation_with_fitted_lines(prep_perc_example$df, prep_perc_example$line_before, prep_perc_example$line_after, "pdist", "normalized pdist", threshold)+theme(legend.position = "right")
hexbin_example_imd <- hexbin_example_imd+theme(legend.position = "right")
hexbin_example_tm <-hexbin_example_tm+theme(legend.position = "none")
hexbin_example_perc <- hexbin_example_perc+theme(legend.position = "none")


# -----------------------------------------------------------------------------
#  MAKE PANEL 1
# -----------------------------------------------------------------------------


# R2 plots
pr2_before <- plot_distribution_across_families(slopes_per_family, "r2_before")+theme(legend.position = "none")+xlab("R² before threshold")
pr2_after <- plot_distribution_across_families(slopes_per_family, "r2_after")+theme(legend.position = "none")+xlab("R² after threshold")

# Slopes plots
p_b_before <- plot_distribution_across_families(slopes_per_family, "before")+theme(legend.position = "none")+xlab("Slope before threshold")
p_b_after <- plot_distribution_across_families(slopes_per_family, "after")+theme(legend.position = "right")+xlab("Slope after threshold")

# Delta
b_delta<- plot_boxplot_across_families(slopes_per_family, "diff")+ylab("Delta slope before vs after")+theme(legend.position = "none")+scale_y_log10()
b_delta


empty_panel <- plot_spacer()
c1 <- p_sat1 + hexbin_example_perc  + plot_annotation(tag_levels = list("A", "B"))
c1b <-  hexbin_example_tm + hexbin_example_imd + plot_annotation(tag_levels = list("C", "D"))
c2 <- pr2_before + p_b_before + b_delta + plot_annotation(tag_levels = list("E", "F", "G"))
c3 <- pr2_after + p_b_after  +empty_panel+ plot_annotation(tag_levels = list("H", "I"))
combined_plot <- (c1 /c1b/ c2 / c3) + plot_annotation(tag_levels = 'A')
ggsave(paste(source_data, '../plots/review/', "SATURATION.png", sep = ''), plot = combined_plot, width = 10, height = 15, dpi = 300, units = 'in')







# # -----------------------------------------------------------------------------
# #  GET EXAMPLES
# # -----------------------------------------------------------------------------
# plot_scatters_distances <- function(input_patristic_df, low_5, plot_name){
#   # plot the low families
#   plot_low_5_list = list()
#   i = 1 
#   for (fam in low_5){

    
#     # compute R
#     input_patristic_df_subset = input_patristic_df[input_patristic_df$family == fam,]

#     # compute R 
#     r = cor(input_patristic_df_subset$dIMD, input_patristic_df_subset$pML, method = "pearson")
#     r2 = r^2
#     plot_sat = saturation_plot_smooth(input_patristic_df_subset, "pML", "dIMD", "ML patristic distance", "IMD distance")
#     plot_sat = plot_sat +ggtitle(paste(fam, "\n R² = ", round(r2, 2), "\n pdist = ", round(mean(input_patristic_df_subset$oneminperc), 2)))
#     # add pdist as title too 
#     # put title in the milddle
#     plot_sat = plot_sat + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
#     #ggsave(paste(paste(plots, "/families_low/", sep = ""),'_',fam,'_IMD_ML.png', sep=''), plot=plot_sat)

  
#     plot_low_5_list[[i]] = plot_sat 

#     # also add tm ml
#     # compute R
#     r = cor(input_patristic_df_subset$dTM, input_patristic_df_subset$pML, method = "pearson")
#     r2 = r^2
#     plot_sat = saturation_plot_smooth(input_patristic_df_subset, "pML", "dTM", "ML patristic distance", "TM distance")+ggtitle(paste("R² = ", round(r2, 2)))
#     plot_sat = plot_sat +ggtitle(paste(fam, "\n R² = ", round(r2, 2), "\n pdist = ", round(mean(input_patristic_df_subset$oneminperc), 2)))
#     # put title in the milddle
#     plot_sat = plot_sat + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
#     plot_low_5_list[[i+1]] = plot_sat

#     # add perc id
#     r = cor(input_patristic_df_subset$oneminperc, input_patristic_df_subset$pML, method = "pearson")
#     r2 = r^2
#     plot_sat = saturation_plot_smooth(input_patristic_df_subset, "pML", "oneminperc", "ML patristic distance", "pdist")+ggtitle(paste("R² = ", round(r2, 2)))
#     plot_sat = plot_sat +ggtitle(paste(fam, "\n R² = ", round(r2, 2), "\n pdist = ", round(mean(input_patristic_df_subset$oneminperc), 2)))
#     # put title in the milddle
#     plot_sat = plot_sat + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
#     plot_low_5_list[[i+2]] = plot_sat
#     i = i + 3

#     # save input patristic subset
#     # save source data

#   }

#   # put the low 5 in a panel with 2 col using patchwork
#   # Combine plots into a grid with one common legend
#   combined_plot <- (plot_low_5_list[[1]]+
#                     plot_low_5_list[[2]] + 
#                     plot_low_5_list[[3]] +
#                     plot_low_5_list[[4]] +
#                     plot_low_5_list[[5]] +
#                     plot_low_5_list[[6]] +
#                     plot_low_5_list[[7]] +
#                     plot_low_5_list[[8]] +
#                     plot_low_5_list[[9]] +
#                     plot_low_5_list[[10]]+
#                     plot_low_5_list[[11]]+
#                     plot_low_5_list[[12]]+
#                     plot_low_5_list[[13]]+
#                     plot_low_5_list[[14]]+
#                     plot_low_5_list[[15]]
#                     ) + 
#     plot_layout(ncol = 3, guides = "collect") & 
#     theme(legend.position = 'right')&
#     plot_annotation(tag_levels = 'A')

#   # save the final image
#   ggsave(paste(source_data, '../plots/review/', plot_name, sep = ''), plot = combined_plot, width = 12, height = 18, dpi = 300, units = 'in')

# }

# # sort the data frame by cor, lowest first
# input_patristic_df$oneminperc = 100 - input_patristic_df$percid

# # Get bottom at 9th decile
# perc_quantile_threshold_10 = quantile(correlations$cor_dIMD_pML, c(0.1)) 
# removed_quantile_10 = correlations[correlations$cor_dIMD_pML > perc_quantile_threshold_10,]
# removed_quantile_10 = removed_quantile_10[order(removed_quantile_10$cor_dIMD_pML),]
# low_5 = removed_quantile_10[1:5,"family"]
# low_5

# # Get middle at 5th decile
# perc_quantile_threshold_50 = quantile(correlations$cor_dIMD_pML, c(0.5, 0.4))
# # get the values in the middle
# quantile_50_elements = correlations[correlations$cor_dIMD_pML > perc_quantile_threshold_50[2] & correlations$cor_dIMD_pML < perc_quantile_threshold_50[1],]
# # now get the middle ones, sort and get the middle
# # order 
# quantile_50_elements = quantile_50_elements[order(quantile_50_elements$cor_dIMD_pML),]
# middle_5th_idx = round(nrow(quantile_50_elements)/2)
# lower_middle_5 = middle_5th_idx - 2
# upper_middle_5 = middle_5th_idx + 2
# mid_5 = quantile_50_elements[lower_middle_5:upper_middle_5,"family"]
# mid_5

# # get 5 random elements near the median 
# # extract the median
# med = median(correlations$cor_dIMD_pML)
# # get the 5 closest values to the median
# correlations$diff = abs(correlations$cor_dIMD_pML - med)
# # sort by diff
# correlations = correlations[order(correlations$diff),]
# mid_5 = correlations[1:5,]$family

# # now extract the 5 families with the most number of sequences
# # get the number of sequences per family
# seqs_per_fam = data.frame()
# for (fam in fl){
#   seqs_per_fam = rbind(seqs_per_fam, data.frame(family = fam, n_seqs = nrow(input_patristic_df[input_patristic_df$family == fam,])))
# }

# # get the 5 families with the most sequences
# top_5_powered = seqs_per_fam[order(seqs_per_fam$n_seqs, decreasing = TRUE),][1:5,"family"]
# top_5_powered

# # get the ones with the highest correlation
# top_5 = correlations[order(correlations$cor_dIMD_pML, decreasing = TRUE),][1:5,"family"]



# #--------------------------------------------------------------------------
# # PLOT THE SCATTERS
# #--------------------------------------------------------------------------
# plot_scatters_distances(input_patristic_df, low_5, 'SATURATION_bottom_9th_decile.png')
# plot_scatters_distances(input_patristic_df, mid_5, 'SATURATION_middle_5th_decile.png')
# plot_scatters_distances(input_patristic_df, top_5_powered, 'SATURATION_top_5_powered.png')
# plot_scatters_distances(input_patristic_df, top_5, 'SATURATION_top_5.png')

# #--------------------------------------------------------------------------
# # SAVE SOURCE DATA
# #--------------------------------------------------------------------------
# input_patristic_df_for_sd = input_patristic_df[,c("family", "seqs", "dIMD", "dTM", "pML", "oneminperc")]

# # save the input patristic distances for low_5 families in a file
# write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% low_5,], paste(source_data, 'saturation_low_5.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")

# # save the input patristic distances for mid_5 families in a file
# write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% mid_5,], paste(source_data, 'saturation_mid_5.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")
# # save the input patristic distances for top_5 families in a file
# write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% top_5,], paste(source_data, 'saturation_top_5.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")
# # save the input patristic distances for top_5 families in a file
# write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% top_5_powered,], paste(source_data, 'saturation_top_5_powered.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")