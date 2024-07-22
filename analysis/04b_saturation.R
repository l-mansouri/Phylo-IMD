library(phangorn)
library(ggplot2)
library(ggExtra)
library(reshape2)
aligners=c('mTMalign')
trimming=c('untrimmed')
source_data = '/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data/'
plots = paste(source_data, '../plots/', sep = '')

fl=read.table(paste(source_data, 'list_of_families_with_all_rep_in_3d', sep = '/'))[,1]

al=aligners[1]
tr=trimming[1]

# read in data frame
input_patristic_df = read.table(paste(source_data,al,'_',tr,'_input_patristic_distances_test.txt', sep=''), header = TRUE)
perc_sim = read.table(paste(source_data, 'percids.txt', sep = ''), header = TRUE)
# create a new column with seqquence 1 and 2 ordered in alphaebtical order
perc_sim$id = pmin(perc_sim$seq1, perc_sim$seq2)
# now append the max 
perc_sim$id = paste(perc_sim$id, pmax(perc_sim$seq1, perc_sim$seq2), sep = "_")
# now append family name
perc_sim$id = paste(perc_sim$id, perc_sim$family, sep = "_")
# subsittue .pdb with .p 
perc_sim$id = gsub(".pdb", ".p", perc_sim$id)
# remove family column


# now do the same for the input_patristic_df
# create seq1 and seq2 columns by splitting the seqs column 
input_patristic_df$seq1 = sapply(strsplit(input_patristic_df$seqs, "_"), "[", 1)
input_patristic_df$seq2 = sapply(strsplit(input_patristic_df$seqs, "_"), "[", 2)
input_patristic_df$id = pmin(input_patristic_df$seq1, input_patristic_df$seq2)
input_patristic_df$id = paste(input_patristic_df$id, pmax(input_patristic_df$seq1, input_patristic_df$seq2), sep = "_")
input_patristic_df$id = paste(input_patristic_df$id, input_patristic_df$family, sep = "_")
# remove seq1 and seq2 columns
input_patristic_df = input_patristic_df[, !(names(input_patristic_df) %in% c("seq1", "seq2"))]

# keep only first family 
# fl = fl[1]
# perc_sim = perc_sim[perc_sim$family == fl,]
# input_patristic_df = input_patristic_df[input_patristic_df$family == fl,]
perc_sim = perc_sim[, !(names(perc_sim) %in% c("family"))]
# now merge the two data frames
input_patristic_df = merge(input_patristic_df, perc_sim, by = "id", all.x = TRUE)

# # Local correction 
# # now per family, apply correction by dividing by the row with max pIMD
# for (fam in fl){
#   input_patristic_df_subset = input_patristic_df[input_patristic_df$family == fam,]
#   max_pIMD = input_patristic_df_subset[which.max(input_patristic_df_subset$pML),]
#   input_patristic_df_subset$pIMD = input_patristic_df_subset$pIMD/max_pIMD$pIMD
#   input_patristic_df_subset$pME = input_patristic_df_subset$pME/max_pIMD$pME
#   input_patristic_df_subset$pTM = input_patristic_df_subset$pTM/max_pIMD$pTM
#   input_patristic_df_subset$pML = input_patristic_df_subset$pML/max_pIMD$pML
#   input_patristic_df_subset$dIMD = input_patristic_df_subset$dIMD/max_pIMD$dIMD
#   input_patristic_df_subset$dME = input_patristic_df_subset$dME/max_pIMD$dME
#   input_patristic_df[input_patristic_df$family == fam,] = input_patristic_df_subset
# }


# for (fam in fl){
#   max_pML = input_patristic_df[input_patristic_df$family == fam,][which.max(input_patristic_df[input_patristic_df$family == fam,]$pML),]
#   max_pML = max_pML$pML
#   input_patristic_df_subset = input_patristic_df[input_patristic_df$family == fam,]
#   input_patristic_df_subset$pML = input_patristic_df_subset$pML/max_pML
#   input_patristic_df_subset$dIMD = input_patristic_df_subset$dIMD/max_pML
#   input_patristic_df[input_patristic_df$family == fam,] = input_patristic_df_subset
# }

# # global correction 
# max_pIMD = input_patristic_df[which.max(input_patristic_df$pML),]
# input_patristic_df$pIMD = input_patristic_df$pIMD/max_pIMD$pIMD


# global correction 
# max_pIMD = input_patristic_df[which.max(input_patristic_df$pIMD),]
# input_patristic_df$pIMD = input_patristic_df$pIMD/max_pIMD$pIMD
# input_patristic_df$pME = input_patristic_df$pME/max_pIMD$pME
# input_patristic_df$pTM = input_patristic_df$pTM/max_pIMD$pTM
# input_patristic_df$pML = input_patristic_df$pML/max_pIMD$pML
# input_patristic_df$dIMD = input_patristic_df$dIMD/max_pIMD$dIMD
# input_patristic_df$dME = input_patristic_df$dME/max_pIMD$dME



saturation_plot <- function(input_patristic_df, distance, patristic, distance_label = "", patristic_label = ""){
  max_x = max(input_patristic_df[,distance])
  max_y = max(input_patristic_df[,patristic])
  max_lim = max(max_x, max_y)+1
  max_lim = 20 
  plot_saturation =ggplot(input_patristic_df, aes_string(x=distance, y=patristic))+
          geom_point(alpha=0)+
          geom_hex(bins=50)+
          xlab(paste(distance_label, ' distance', sep = " "))+
          ylab(paste(patristic_label, '', sep = " "))+
          theme_light()+
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
          geom_point(alpha=0.1)+
          geom_smooth(method = "lm")+
          xlab(paste(distance_label, ' distance', sep = " "))+
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


#filter out entries above 6 patristic
#input_patristic_df = input_patristic_df[input_patristic_df$pML < 7,]
#input_patristic_df = input_patristic_df[input_patristic_df$dIMD < 3,]

head(input_patristic_df)
# remove the rows with all 0s

saturation_plot(input_patristic_df, "dME", "pME", "ME", "ME")
saturation_plot(input_patristic_df, "dIMD", "pIMD", "IMD", "IMD")
saturation_plot(input_patristic_df, "dTM", "pTM", "TM", "TM")
saturation_plot(input_patristic_df, "dIMD", "pML", "IMD", "ML")
saturation_plot(input_patristic_df, "dME", "pML", "ME", "ML")

saturation_plot(input_patristic_df, "dIMD", "percid", "IMD", "percid")
saturation_plot(input_patristic_df, "dTM", "percid", "TM", "percid")

# subset the plot for one family 

for (fam in fl){
  input_patristic_df_subset = input_patristic_df[input_patristic_df$family == fam,]
  plot_sat = saturation_plot_smooth(input_patristic_df_subset, "dIMD", "pML", "IMD", "ML")
  # plot_me = saturation_plot_smooth(input_patristic_df_subset, "dME", "pML", "ME", "ML")
  # plot_tm = saturation_plot_smooth(input_patristic_df_subset, "dTM", "pML", "TM", "ML")
  # plot_imd_tm = saturation_plot_smooth(input_patristic_df_subset, "dIMD", "pTM", "IMD", "TM")
  
  # save the plot 
  ggsave(paste(paste(plots, "/families/", sep = ""),al,'_',tr,'_',fam,'_IMD_ML_saturation_smooth_corrected.png', sep=''), plot=plot_sat)
  # ggsave(paste(paste(plots, "/families/", sep = ""),al,'_',tr,'_',fam,'_ME_ML_saturation_smooth_corrected.png', sep=''), plot=plot_me)
  # ggsave(paste(paste(plots, "/families/", sep = ""),al,'_',tr,'_',fam,'_TM_ML_saturation_smooth_corrected.png', sep=''), plot=plot_tm)
  # ggsave(paste(paste(plots, "/families/", sep = ""),al,'_',tr,'_',fam,'_IMD_TM_saturation_smooth_corrected.png', sep=''), plot=plot_imd_tm)
}




# Compute the correlation between the distances and the patristic distances per family 

correlations = data.frame()
for (fam in fl){
  input_patristic_df_subset = input_patristic_df[input_patristic_df$family == fam,]
  cor_dTM_pTM = cor(input_patristic_df_subset$dTM, input_patristic_df_subset$pTM)
  cor_dIMD_pIMD = cor(input_patristic_df_subset$dIMD, input_patristic_df_subset$pIMD)
  cor_dME_pME = cor(input_patristic_df_subset$dME, input_patristic_df_subset$pME)
  cor_dIMD_pML = cor(input_patristic_df_subset$dIMD, input_patristic_df_subset$pML)
  correlations = rbind(correlations, data.frame(family = fam, cor_dTM_pTM = cor_dTM_pTM, cor_dIMD_pIMD = cor_dIMD_pIMD, cor_dME_pME = cor_dME_pME, cor_dIMD_pML = cor_dIMD_pML))
}

correlations
# remove family column and square the rest
correlations = correlations[, !(names(correlations) %in% c("family"))]
correlations = correlations^2


# plot distribution cor_dIMD_pML
ggplot(mtm, aes(x = COR3dml, fill = "grey"))+
  geom_density()+
  xlab('R (input distances IMD, patristic distances ML)')+
  ylab('Density')+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black", size = 10),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),    # Hide axis ticks
        panel.grid = element_blank(), 
        plot.background = element_rect(fill = "white"),
        axis.ticks.y = element_line(color = "black")) +
  expand_limits(y=0)

# square correlations


correlations


# plot the correlations - boxplot
correlations_melt = melt(correlations, id.vars = "family")
ggplot(correlations_melt, aes(x = variable, y = value, fill = variable))+
  geom_boxplot(outlier.size=0.1, lwd=0.6, linetype = "solid" )+
  xlab('')+
  ylab('R² (input distances, patristic distances)')+
  theme_minimal()+
  theme(
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_blank(),    # Hide axis text
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.ticks.x = element_blank(),    # Hide axis ticks
        panel.grid = element_blank(), 
        plot.background = element_rect(fill = "white"),
        axis.ticks.y = element_line(color = "black")) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1), linetype = "solid", color = "black")+
  expand_limits(y=0)

median_seq = round(median(correlations$cor_dME_pME),3)
median_imd = round(median(correlations$cor_dIMD_pIMD),3)
median_tm  = round(median(correlations$cor_dTM_pTM),3)

median_seq
median_imd
median_tm


mtm
