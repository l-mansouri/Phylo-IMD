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

# find outliers of using quantiles
qdimd = quantile(input_patristic_df$dIMD, 0.995)
qml = quantile(input_patristic_df$pML, 0.995)
qtm = quantile(input_patristic_df$dTM, 0.995)


# Calculate the slope before and after the 0.25 threshold om the ML patristic distance
columns_filter = c("id", "dIMD", "dTM", "pML", "family", "percid", "seq1", "seq2")
df = input_patristic_df[,columns_filter]
df$pdist = 100-df$percid

# -----------------------------------------------------------------------------
#       IDENTIFY TWILIGHT ZONE THRESHOLD
# -----------------------------------------------------------------------------

# get the threshold by finding the point of the line with perc_identity =70% 
df_percid_ml = df[,c("pdist", "pML")]
df_percid_ml = df_percid_ml[df_percid_ml$pML < qml,]
ml_perc <- saturation_plot(df_percid_ml,  "pML","pdist", "ML patristic distance", "pdist")+geom_smooth()
fitted_values <- ggplot_build(ml_perc)$data[[5]]
fitted_values = fitted_values[fitted_values$y <= 70,]
max_y = max(fitted_values$y)
threshold = fitted_values[fitted_values$y == max_y,]$x
threshold

# -----------------------------------------------------------------------------
#       TABLE TWILIGHT ZONE
# -----------------------------------------------------------------------------


prep_df_for_plotting <- function(df, y, quantile_threshold, threshold, apply_filter = TRUE, normalization_factor = NA){
  # IMD vs ML patristic distance
  df_imd_ml = df[,c(y, "pML", "family", "seq1", "seq2")]
  if(apply_filter){
    print("FILTERING")
    df_imd_ml = df_imd_ml[df_imd_ml[[y]] < quantile_threshold & df_imd_ml$pML < qml,]
  }
  var_norm = paste(y, "_norm", sep = "")
  if(!is.na(normalization_factor)){
    df_imd_ml[[var_norm]] = df_imd_ml[[y]]*100/normalization_factor
  }else{
    df_imd_ml[[var_norm]] = df_imd_ml[[y]]*100/median(df_imd_ml[[y]])
  }
  df_imd_ml_before = df_imd_ml[df_imd_ml$pML < threshold,]
  df_imd_ml_after = df_imd_ml[df_imd_ml$pML >= threshold,]
  # Check if there are enough points to fit a line
  if(nrow(df_imd_ml_before) < 2 | nrow(df_imd_ml_after) < 2){
    return( list(df = df_imd_ml, line_before = NA, line_after = NA, slope_before = NA, slope_after = NA, r2_before = NA, r2_after = NA))
  }
  line_imd_ml_before = lm(as.formula(paste(var_norm, "~", "pML")), data = df_imd_ml_before)
  slope_imd_ml_before <- coef(line_imd_ml_before)[2]
  line_imd_ml_after = lm(as.formula(paste(var_norm, "~", "pML")), data = df_imd_ml_after)
  slope_imd_ml_after <- coef(line_imd_ml_after)[2]
  # round
  slope_imd_ml_before = round(slope_imd_ml_before, 2)
  slope_imd_ml_after = round(slope_imd_ml_after, 2)
  # calculate R²
  r_imd_ml_before = cor(df_imd_ml_before[[var_norm]], df_imd_ml_before$pML, method = "pearson")
  r2_imd_ml_before = round(r_imd_ml_before^2,2)
  r_imd_ml_after = cor(df_imd_ml_after[[var_norm]], df_imd_ml_after$pML, method = "pearson")
  r2_imd_ml_after = round(r_imd_ml_after^2,2)

  # return the data frame, line_imd_ml_before, line_imd_ml_after, slope_imd_ml_before, slope_imd_ml_after, r2_imd_ml_before, r2_imd_ml_after
  return( list(df = df_imd_ml, line_before = line_imd_ml_before, line_after = line_imd_ml_after, slope_before = slope_imd_ml_before, slope_after = slope_imd_ml_after, r2_before = r2_imd_ml_before, r2_after = r2_imd_ml_after))

}


imd_prep  <- prep_df_for_plotting(df, "dIMD", qdimd, threshold)
df_imd_ml <- imd_prep$df
tm_prep  <- prep_df_for_plotting(df, "dTM", qtm, threshold)
df_tm_ml <- tm_prep$df
perc_prep  <- prep_df_for_plotting(df, "pdist", Inf, threshold)
df_percid_ml <- perc_prep$df




slopes = data.frame(method = c("IMD ", "TM ", "100-%id"),
                   before = c(imd_prep$slope_before, tm_prep$slope_before, perc_prep$slope_before),
                   after = c(imd_prep$slope_after, tm_prep$slope_after, perc_prep$slope_after),
                    r2_before = c(imd_prep$r2_before, tm_prep$r2_before, perc_prep$r2_before),
                    r2_after = c(imd_prep$r2_after, tm_prep$r2_after, perc_prep$r2_after))
slopes$ratio = round(slopes$before/slopes$after, 2)
write.table(slopes, paste(source_data, '../tables/slopes_twilight_zone.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")
slopes
# make the slopes table without filtering
imd_prep_unfiltered  <- prep_df_for_plotting(df, "dIMD", qdimd, threshold, apply_filter = FALSE)
tm_prep_unfiltered  <- prep_df_for_plotting(df, "dTM", qtm, threshold, apply_filter = FALSE)
perc_prep_unfiltered  <- prep_df_for_plotting(df, "pdist", Inf, threshold, apply_filter = FALSE)
slopes_unfiltered = data.frame(method = c("IMD ", "TM ", "100-%id"),
                   before = c(imd_prep_unfiltered$slope_before, tm_prep_unfiltered$slope_before, perc_prep_unfiltered$slope_before),
                   after = c(imd_prep_unfiltered$slope_after, tm_prep_unfiltered$slope_after, perc_prep_unfiltered$slope_after),
                    r2_before = c(imd_prep_unfiltered$r2_before, tm_prep_unfiltered$r2_before, perc_prep_unfiltered$r2_before),
                    r2_after = c(imd_prep_unfiltered$r2_after, tm_prep_unfiltered$r2_after, perc_prep_unfiltered$r2_after))
slopes_unfiltered$ratio = round(slopes_unfiltered$before/slopes_unfiltered$after, 2)
write.table(slopes_unfiltered, paste(source_data, '../tables/slopes_twilight_zone_unfiltered.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")



# for each family, split at threshold and caluculate r2 and slope before and after
# normalize by the median 

filter = TRUE
if(filter){
  norm_imd = median(df_imd_ml$dIMD)
  norm_tm = median(df_tm_ml$dTM)
  norm_perc = median(df_percid_ml$pdist)
}else{
  norm_imd = median(df$dIMD)
  norm_tm = median(df$dTM)
  norm_perc = median(df$pdist)
}

slopes_per_family = data.frame()
threshold = median(df$pML)
print(threshold)
for (fam in fl){
  df_subset = df[df$family == fam,]
  imd_subset_prep = prep_df_for_plotting(df_subset, "dIMD", qdimd, threshold, apply_filter = filter, normalization_factor = norm_imd)

  tm_subset_prep = prep_df_for_plotting(df_subset, "dTM", qtm, threshold, apply_filter = filter, normalization_factor = norm_tm)

  perc_subset_prep = prep_df_for_plotting(df_subset, "pdist", Inf, threshold, apply_filter = filter, normalization_factor = norm_perc)

  slopes = data.frame(method = c("IMD ", "TM ", "pdist"),
                   before = c(imd_subset_prep$slope_before, tm_subset_prep$slope_before, perc_subset_prep$slope_before),
                   after = c(imd_subset_prep$slope_after, tm_subset_prep$slope_after, perc_subset_prep$slope_after),
                    r2_before = c(imd_subset_prep$r2_before, tm_subset_prep$r2_before, perc_subset_prep$r2_before),
                    r2_after = c(imd_subset_prep$r2_after, tm_subset_prep$r2_after, perc_subset_prep$r2_after))

  slopes$ratio = round(slopes$before/slopes$after, 2)
  slopes$family = fam
  slopes_per_family = rbind(slopes_per_family, slopes)
}

# Check how many families i have removed from the analysis
families_excluded <- (unique(slopes_per_family$family[is.na(slopes_per_family$before)]))


# Plot distribution of R² values before threshold, colored by method
# remove NAs
slopes_per_family = slopes_per_family[!is.na(slopes_per_family$before),]

plot_distribution_across_families <- function(slopes_per_family, y){
  correlations_melted = melt(slopes_per_family[,c(y, "method")])
  # plot the distribution of R² values before
  correlations_melted$method = factor(correlations_melted$method, levels = c("IMD ", "TM ", "pdist"))
  slopes_per_family$y = slopes_per_family[[y]]
  medians = aggregate(y ~ method, data = slopes_per_family, FUN = median)
  p <- ggplot(correlations_melted, aes(x = value, fill = method, color = method)) +
              geom_density(alpha = 0.5) +
              theme_light() +
              xlab(y) + 
              ylab("Density") +
              theme(legend.title=element_blank())+ theme(legend.position = "right")+
              geom_vline(data = medians, aes(xintercept = y, color = method), linetype = "dashed", linewidth = 1)
  return(p)
}

plot_boxplot_across_families <- function(slopes_per_family, y){
  correlations_melted = melt(slopes_per_family[,c(y, "method")])
  # plot the distribution of R² values before
  correlations_melted$method = factor(correlations_melted$method, levels = c("IMD ", "TM ", "pdist"))
  slopes_per_family$y = slopes_per_family[[y]]
  p <- ggplot(correlations_melted, aes(x = method, y = value, fill = method, color = method)) +
              geom_boxplot(alpha = 0.5) +
              theme_light() +
              xlab("") + 
              ylab(paste(y, " (log10)")) +
              # log scale
              scale_y_log10() +
              theme(legend.title=element_blank())+ theme(legend.position = "right")
  return(p)
}


pr2_before <- plot_distribution_across_families(slopes_per_family, "r2_before")
pr2_after <- plot_distribution_across_families(slopes_per_family, "r2_after")
p_b_before <- plot_distribution_across_families(slopes_per_family, "before")
b_ratio <- plot_boxplot_across_families(slopes_per_family, "ratio")

c1 <- pr2_before + pr2_after + plot_annotation(tag_levels = list("A", "B"))
c2 <- p_b_before + b_ratio + plot_annotation(tag_levels = list("C", "D"))
combined_plot <- (c1 / c2) + plot_annotation(tag_levels = 'A')
# add title 
# if filter is TRUE, filtering criteria is "filtered", otherwise "unfiltered"
filtering_criteria = ifelse(filter, "filtered", "unfiltered")
combined_plot = combined_plot + ggtitle(filtering_criteria)
name_figure = paste("panel_SAT_", filtering_criteria, ".png", sep = "")
print(name_figure)
ggsave(paste(source_data, '../plots/review/', name_figure, sep = ''), plot = combined_plot, width = 10, height = 8, dpi = 300, units = 'in')


# -----------------------------------------------------------------------------
#       PLOT PANEL
# -----------------------------------------------------------------------------

maxfill_pdist <- max(hexbin(df_percid_ml[["pML"]], df_percid_ml[["pdist_norm"]], xbins = 50)@count)
maxfill_imd <- max(hexbin(df_imd_ml[["pML"]], df_imd_ml[["dIMD_norm"]], xbins = 50)@count)
maxfill_tm <- max(hexbin(df_tm_ml[["pML"]], df_tm_ml[["dTM_norm"]], xbins = 50)@count)
maxfillval = max(maxfill_pdist, maxfill_imd, maxfill_tm)
maxfillval



# -----------------------------------------------------------------------------
#       PERC ID
# -----------------------------------------------------------------------------


plot_saturation_with_fitted_lines <- function(df_imd_ml,line_imd_ml_before, line_imd_ml_after, y, ylabel,  threshold){
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


# PLOT
pml_dimd <- plot_saturation_with_fitted_lines(df_imd_ml,line_imd_ml_before, line_imd_ml_after, "dIMD_norm","IMD distance", threshold)
pml_dtm  <- plot_saturation_with_fitted_lines(df_tm_ml,line_tm_ml_before, line_tm_ml_after, "dTM_norm", "TM distance", threshold)
ml_perc  <- plot_saturation_with_fitted_lines(df_percid_ml,line_percid_ml_before, line_percid_ml_after, "pdist_norm", "pdist", threshold)
ml_perc
# COMBINE PLOTS
c1 <- p_sat1 + ml_perc + theme(legend.position = "none") + plot_annotation(tag_levels = 'A')
c2 <- pml_dtm + pml_dimd + theme(legend.position = "right") + plot_annotation(tag_levels = list("C", "D"))
combined_plot <- (c1 / c2) + plot_annotation(tag_levels = 'A')
ggsave(paste(source_data, '../plots/review/', 'SATURATION_plots.png', sep = ''), plot = combined_plot, width = 13, height = 11, dpi = 300, units = 'in')



# -----------------------------------------------------------------------------
#  COMPUTE CORRELATIONS
# -----------------------------------------------------------------------------
correlations = data.frame()
for (fam in fl){
  df_imd_ml_subset = df_imd_ml[df_imd_ml$family == fam,]
  cor_dIMD_pML = cor(df_imd_ml_subset$dIMD, df_imd_ml_subset$pML, method = "pearson")
  r2_imd_ml = cor_dIMD_pML^2

  df_tm_ml_subset = df_tm_ml[df_tm_ml$family == fam,]
  cor_dTM_pML = cor(df_tm_ml_subset$dTM, df_tm_ml_subset$pML, method = "pearson")
  r2_tm_ml = cor_dTM_pML^2

  # perc id make it 100 - perc id
  df_percid_ml_subset = df_percid_ml[df_percid_ml$family == fam,]
  cor_percid_pML = cor(df_percid_ml_subset$pdist, df_percid_ml_subset$pML, method = "pearson")
  r2_percid_ml = cor_percid_pML^2

  correlations = rbind(correlations, data.frame(family = fam, cor_dIMD_pML = cor_dIMD_pML, r2_imd_ml = r2_imd_ml, cor_dTM_pML = cor_dTM_pML, r2_tm_ml = r2_tm_ml, cor_percid_pML = cor_percid_pML, r2_percid_ml = r2_percid_ml))
}





# plot density of the r2 values
correlations_melted = melt(correlations[,c("r2_imd_ml", "r2_tm_ml", "r2_percid_ml")])
colnames(correlations_melted) = c("method", "r2")
correlations_melted$method = factor(correlations_melted$method, levels = c("r2_imd_ml", "r2_tm_ml", "r2_percid_ml"))
correlations_melted$method = factor(correlations_melted$method, labels = c("IMD vs pat ML", "TM vs pat ML", "pdist vs pat ML"))
# calc medians
medians = aggregate(r2 ~ method, data = correlations_melted, FUN = median)
r2_plot <- ggplot(correlations_melted, aes(x = r2, fill = method, color =method)) +
            geom_density(alpha = 0.5) +
            theme_light() +
            xlab("R²") + 
            ylab("Density") +
            geom_vline(data = medians, aes(xintercept = r2, color = method), linetype = "dashed", linewidth = 1)+
            theme(legend.title=element_blank())+ theme(legend.position = "right")


# -----------------------------------------------------------------------------
#  PLOT UNFILTERED
# -----------------------------------------------------------------------------
# Show filtering lines
df_imd_ml = df[,c("dIMD", "pML", "family", "seq1", "seq2")]
df_tm_ml = df[,c("dTM", "pML", "family", "seq1", "seq2")]
df_percid_ml = df[,c("pdist", "pML", "family", "seq1", "seq2")]
df_percid_ml$percid = (100-df_percid_ml$percid)

maxfill_pdist <- max(hexbin(df_percid_ml[["pML"]], df_percid_ml[["pdist"]], xbins = 50)@count)
maxfill_imd <- max(hexbin(df_imd_ml[["pML"]], df_imd_ml[["dIMD"]], xbins = 50)@count)
maxfill_tm <- max(hexbin(df_tm_ml[["pML"]], df_tm_ml[["dTM"]], xbins = 50)@count)
maxfillval = max(maxfill_pdist, maxfill_imd, maxfill_tm)

pml_dimd_show_filtering <- saturation_plot(df_imd_ml, "pML", "dIMD", "ML patristic distance", "IMD distance", maxfillval)+
  geom_vline(xintercept = qml, color = "red", linetype = "dashed", linewidth = 0.7)+
  geom_hline(yintercept = qdimd, color = "red", linetype = "dashed", linewidth = 0.7)

pml_dtm_show_filtering <- saturation_plot(df_tm_ml, "pML", "dTM", "ML patristic distance", "TM distance", maxfillval)+
  geom_vline(xintercept = qml, color = "red", linetype = "dashed", linewidth = 0.7)+
  geom_hline(yintercept = qtm, color = "red", linetype = "dashed", linewidth = 0.7)+theme(legend.position = "none")

pml_percid_show_filtering <- saturation_plot(df_percid_ml, "pML", "pdist", "ML patristic distance", "pdist", maxfillval)+
  geom_vline(xintercept = qml, color = "red", linetype = "dashed", linewidth = 0.7)+ theme(legend.position = "none")

c1 <-  pml_percid_show_filtering + pml_dimd_show_filtering   + plot_annotation(tag_levels = list("A", "B"))
c2 <-  pml_dtm_show_filtering + r2_plot  + plot_annotation(tag_levels = list("C", "D"))
combined_plot <- (c1 / c2) + plot_annotation(tag_levels = 'A')+theme(legend.position = "rigth")

#save
ggsave(paste(source_data, '../plots/review/', 'supp_saturation.png', sep = ''), width = 10, height = 8, dpi = 300, units = 'in')


# -----------------------------------------------------------------------------
#  GET EXAMPLES
# -----------------------------------------------------------------------------
plot_scatters_distances <- function(input_patristic_df, low_5, plot_name){
  # plot the low families
  plot_low_5_list = list()
  i = 1 
  for (fam in low_5){

    
    # compute R
    input_patristic_df_subset = input_patristic_df[input_patristic_df$family == fam,]

    # compute R 
    r = cor(input_patristic_df_subset$dIMD, input_patristic_df_subset$pML, method = "pearson")
    r2 = r^2
    plot_sat = saturation_plot_smooth(input_patristic_df_subset, "pML", "dIMD", "ML patristic distance", "IMD distance")
    plot_sat = plot_sat +ggtitle(paste(fam, "\n R² = ", round(r2, 2), "\n pdist = ", round(mean(input_patristic_df_subset$oneminperc), 2)))
    # add pdist as title too 
    # put title in the milddle
    plot_sat = plot_sat + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
    #ggsave(paste(paste(plots, "/families_low/", sep = ""),'_',fam,'_IMD_ML.png', sep=''), plot=plot_sat)

  
    plot_low_5_list[[i]] = plot_sat 

    # also add tm ml
    # compute R
    r = cor(input_patristic_df_subset$dTM, input_patristic_df_subset$pML, method = "pearson")
    r2 = r^2
    plot_sat = saturation_plot_smooth(input_patristic_df_subset, "pML", "dTM", "ML patristic distance", "TM distance")+ggtitle(paste("R² = ", round(r2, 2)))
    plot_sat = plot_sat +ggtitle(paste(fam, "\n R² = ", round(r2, 2), "\n pdist = ", round(mean(input_patristic_df_subset$oneminperc), 2)))
    # put title in the milddle
    plot_sat = plot_sat + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
    plot_low_5_list[[i+1]] = plot_sat

    # add perc id
    r = cor(input_patristic_df_subset$oneminperc, input_patristic_df_subset$pML, method = "pearson")
    r2 = r^2
    plot_sat = saturation_plot_smooth(input_patristic_df_subset, "pML", "oneminperc", "ML patristic distance", "pdist")+ggtitle(paste("R² = ", round(r2, 2)))
    plot_sat = plot_sat +ggtitle(paste(fam, "\n R² = ", round(r2, 2), "\n pdist = ", round(mean(input_patristic_df_subset$oneminperc), 2)))
    # put title in the milddle
    plot_sat = plot_sat + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
    plot_low_5_list[[i+2]] = plot_sat
    i = i + 3

    # save input patristic subset
    # save source data

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
                    plot_low_5_list[[15]]
                    ) + 
    plot_layout(ncol = 3, guides = "collect") & 
    theme(legend.position = 'right')&
    plot_annotation(tag_levels = 'A')

  # save the final image
  ggsave(paste(source_data, '../plots/review/', plot_name, sep = ''), plot = combined_plot, width = 12, height = 18, dpi = 300, units = 'in')

}

# sort the data frame by cor, lowest first
input_patristic_df$oneminperc = 100 - input_patristic_df$percid

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
plot_scatters_distances(input_patristic_df, low_5, 'SATURATION_bottom_9th_decile.png')
plot_scatters_distances(input_patristic_df, mid_5, 'SATURATION_middle_5th_decile.png')
plot_scatters_distances(input_patristic_df, top_5_powered, 'SATURATION_top_5_powered.png')
plot_scatters_distances(input_patristic_df, top_5, 'SATURATION_top_5.png')



#--------------------------------------------------------------------------
# SAVE SOURCE DATA
#--------------------------------------------------------------------------
input_patristic_df_for_sd = input_patristic_df[,c("family", "seqs", "dIMD", "dTM", "pML", "oneminperc")]

# save the input patristic distances for low_5 families in a file
write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% low_5,], paste(source_data, 'saturation_low_5.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")

# save the input patristic distances for mid_5 families in a file
write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% mid_5,], paste(source_data, 'saturation_mid_5.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")
# save the input patristic distances for top_5 families in a file
write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% top_5,], paste(source_data, 'saturation_top_5.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")
# save the input patristic distances for top_5 families in a file
write.table(input_patristic_df_for_sd[input_patristic_df_for_sd$family %in% top_5_powered,], paste(source_data, 'saturation_top_5_powered.csv', sep = ''), row.names = FALSE, quote = FALSE, sep = ",")