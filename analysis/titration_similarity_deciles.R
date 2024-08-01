library(dplyr)
library(ggplot2)
library(dplyr)

source_data = '/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data/'

full_families=read.table(paste(source_data, 'list_of_families_with_all_rep_in_3d', sep = ""))[,1]
titration_families = read.table(paste(source_data, 'list_of_families_for_titration', sep = ''), header = F)[,1]

# load perc id
perc_sim = read.table(paste(source_data, 'percids.txt', sep = ''), header = TRUE)
perc_sim = perc_sim %>% group_by(family) %>% summarise(mean_perc_id = mean(percid))
perc_sim$titration = ifelse(perc_sim$family %in% titration_families, 'titration', 'full')

# ----------------------------------------------------------------------------------------------
# Here check how similar is the distribution of perc id between titration and full families
# ----------------------------------------------------------------------------------------------

# plot highligting titration families
p = ggplot(perc_sim, aes(x = mean_perc_id, fill = titration)) +
        geom_density(alpha = 0.5, colour = NA) +
        theme_minimal() +
        theme(legend.position = 'top') +
        labs(x = 'Mean perc id', y = 'Density') +
        scale_fill_manual(values = c('full' = 'dark grey', 'titration' = 'coral'))+
        theme(axis.text.x = element_text(color = 'black', size = 15),
              axis.text.y = element_text(color = 'black', size = 15),
              axis.title.x = element_text(color = 'black', size = 20),
              axis.title.y = element_text(color = 'black', size = 20),
              legend.text = element_text(color = 'black', size = 15),
              legend.title = element_text(color = 'black', size = 20))+
        expand_limits(x = c(0, 100))+
        xlab("Average sequence similarity (%)")+
        ylab("Density")+
        guides(fill = guide_legend(title = 'Dataset'))



titration_perc = perc_sim[perc_sim$family %in% titration_families,]$mean_perc_id
full_perc = perc_sim[perc_sim$family %in% full_families,]$mean_perc_id
# comparet the two distributions with statitsical test
wilcox.test(titration_perc, full_perc, paired = FALSE)

# ----------------------------------------------------------------------------------------------
# SEPARATE FAMILIES IN DECILES 
# ----------------------------------------------------------------------------------------------

perc_id_full = perc_sim[perc_sim$family %in% titration_families,]
# get deciles of perc id from titration families
breaks <- seq(0, 100, 5)
# do the breaks as deciles
breaks = quantile(perc_id_full$mean_perc_id, probs = seq(0, 1, 0.2))
breaks
perc_id_full$decile = cut(perc_id_full$mean_perc_id, breaks = breaks, include.lowest = TRUE)
head(perc_id_full)
# plot how many per family are in each decile
p = ggplot(perc_id_full, aes(x=decile, fill = titration)) + geom_bar(position = 'dodge') + theme_minimal() + theme(legend.position = 'top') +
 labs(x = 'Decile', y = 'Number of families') + scale_fill_manual(values = c('full' = 'grey', 'titration' = 'coral'))

# Now extract families for each decile
deciles = split(perc_id_full, perc_id_full$decile)
# iterate over deciles and extract families
for (i in 1:length(deciles)){
  decile = deciles[[i]]
  print(paste('Decile:', i))
  decile_name = names(deciles)[i]
  print(decile_name)
  # substitute , with - in decile name
  decile_name = gsub(',', '-', decile_name)
  # store the list of families in a file
  # creae a directory to store the deciles
  if (decile$family %>% length > 0){
    dir.create(paste(source_data, 'titration_deciles/', sep = ''), showWarnings = FALSE)
    write.table(decile$family, paste(source_data, 'titration_deciles/deciles', decile_name, '.txt', sep = ''), row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}


# NOW go and reproduce tables on deciles 
# Table with: decile, number of families in it, number of reference branches, number of columns for 75%

# extract number of total reference branches per family
path_reference_numbers ="/home/luisasantus/Desktop/crg_cluster/newphylo/NF_TMalign/mTMalign/titration_every_5/ref_br_ME+3d+ML"
# read in file and per family extract the number of reference branches (forelast column)
# iterate each file in folder, fam name is in filename fam_ 
files = list.files(path_reference_numbers, pattern = 'PF*', full.names = F)
# iterate through files, extract family name (P*_*) and number of reference branches (forelast column) and store in dataframe
reference_branches = data.frame()
for (file in files){
  fam = strsplit(file, split = '_')[[1]][1]
  # read file
  ref = read.table(paste(path_reference_numbers, file, sep = '/'), header = F)
  # extract number of reference branches
  n_ref = ref[1,ncol(ref)-1]
  n_tot = ref[1,ncol(ref)]
  # store in dataframe
  reference_branches = rbind(reference_branches, data.frame(family = fam, n_ref = n_ref, n_tot = n_tot))
}
reference_branches
perc_sim_titration = perc_sim[perc_sim$family %in% titration_families,]
perc_sim_titration

# decile_name, n_families, n_ref_branches, n_columns
summary_table = data.frame()
for (i in 1:length(deciles)){
  decile = deciles[[i]]
  decile_name = names(deciles)[i]
  if (decile$family %>% length > 0){
    decile_families = decile$family
    n_families = length(decile_families)
    decile_ref_branches = reference_branches[reference_branches$family %in% decile_families,]
    n_ref_branches = sum(decile_ref_branches$n_ref)
    n_tot = sum(decile_ref_branches$n_tot)
    summary_table = rbind(summary_table, data.frame(decile = decile_name, n_families = n_families, n_ref_branches = n_ref_branches, n_tot = n_tot))
  }
}

summary_table
# change comman in decile name to -
summary_table$decile = gsub(',', '-', summary_table$decile)


files = list.files(source_data, pattern = 'deciles.*.txt', full.names = T)
thresholds = c(0.6, 0.75, 0.9)
colnames = c("avg1d", "avg3d", "avgML")


plot_titration <- function(df, labels, palette =  c("#f8766d", "#00ba39", "#619cff")){
  
  colnames(df)=c('ncol','avg1d', 'sd1d', 'avg3d', 'sd3d', 'avgML', 'sdML', 'x')
  dfml=data.frame(n_o_col=c(df$ncol, df$ncol, df$ncol),
                  avg_RF=c(df$avg1d, df$avg3d, df$avgML),
                  sd_RF=c(df$sd1d, df$sd3d, df$sdML),
                  type=factor(rep(labels, each=as.numeric(length(df[,1]))), levels=labels)
  )
  print(nrow(dfml))
  my_palette <- palette

  # Plot using custom palette
  pexp=ggplot(dfml, aes(x=dfml[,1], y=dfml[,2], color=type), alpha=0.2)+
    geom_point()+
    geom_line()+
    theme_minimal()+
    xlab('Number of columns included')+
    ylab('Fraction of reference branches present')+
    geom_ribbon(aes(ymin = dfml[,2] - dfml[,3], ymax = dfml[,2] + dfml[,3], fill = type), color = NA, alpha = 0.2)+
    scale_color_manual(values = my_palette)+
    scale_fill_manual(values = my_palette)+
    xlim(c(0,201))+
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), limits = c(0,1.06))
    # white background
  pexp=pexp+theme(panel.background = element_rect(fill = 'white', colour = 'white'))
  pexp=pexp+geom_hline(yintercept=0.75, linetype="dashed", color = "grey", linewidth=1)

  return(pexp)  
}




# iterate through files, extract file name, save plot with that name
df_thresholds_deciles = data.frame()
plot_list = list()
for (file in files){
  # extract file name
  file_name = gsub('.*/', '', file)
  file_name = gsub('.txt', '', file_name)
  decile = gsub('deciles', '', file_name)
  print(decile)
  decile = strsplit(decile, '_')[[1]][1]
  # read file
  df = read.table(file)
  colnames(df)=c('ncol','avg1d', 'sd1d', 'avg3d', 'sd3d', 'avgML', 'sdML', 'x')
  # plot

  # for each threshold and each column name, get the ncol value and store in df with the corresponding threshold as column and column name as row
  df_thresholds = data.frame(matrix(ncol = length(thresholds), nrow = length(colnames)))
  colnames(df_thresholds) = thresholds
  rownames(df_thresholds) = colnames
  for (threshold in thresholds){
      for (colname in colnames){
          diff =abs(df[colname] - threshold)
          df_thresholds[colname,as.character(threshold)] = df[which.min(diff[[colname]]),"ncol"]
      }
  }

  # rename rows to Seq-ME, IMD-ME, Seq-ML
  rownames(df_thresholds) = c("ME", "IMD", "ML")
  # rename columns to percentage
  colnames(df_thresholds) = c("60%", "75%", "90%")
  # reorder rows
  df_thresholds = df_thresholds[c("ME", "ML","IMD"),]
  df_thresholds$method = rownames(df_thresholds)
  # make methof first column
  df_thresholds = df_thresholds[,c(4,1,2,3)]
  print(df_thresholds)
  p = plot_titration(df, c("ME", "IMD", "ML"))
  # save plot
  # as title add the name of decile and the number of families
  nfam = unique(summary_table[summary_table$decile == decile,]$n_families)
  print(nfam)
  
  p = p + ggtitle(paste('Quantile:', decile, '% \n number of datasets=', nfam))
  # make title centered
  p = p + theme(plot.title = element_text(hjust = 0.5))
  plot_list[[length(plot_list)+1]] = p

  # add decile to the table
  df_thresholds$decile = decile
  df_thresholds_deciles = rbind(df_thresholds_deciles, df_thresholds)
  #ggsave(paste(source_data, '../plots/review/', file_name, '.png', sep = ''), plot = p, width = 7, height = 5.7, dpi = 300, units = 'in')
}

library(patchwork)
plot_list <- c(plot_list[length(plot_list)], plot_list[-length(plot_list)])

add_label_to_plot <- function(plot, label) {
  plot + 
    annotate("text", 
             x = -Inf, y = Inf, 
             label = label, 
             hjust = -0.1, vjust = 1.1, 
             size = 6, fontface = "bold", color = "black") +
    theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 60))
}
# Define labels
labels <- LETTERS[1:length(plot_list)]

# Apply the label function to each plot
plot_list_with_labels <- mapply(add_label_to_plot, plot_list, labels, SIMPLIFY = FALSE)

# Combine plots into a grid
combined_plot <- do.call(grid.arrange, c(plot_list_with_labels, ncol = 3))
# Save the final image
ggsave(paste(source_data, '../plots/review/', 'deciles_plots.png', sep = ''), plot = combined_plot, width = 20, height = 12, dpi = 300, units = 'in')




# -------------------------------------
#       PREPARE SUMMARY TABLE 
# -------------------------------------

summary_table

df_thresholds_deciles


# merge the two tables
summary_table = merge(summary_table, df_thresholds_deciles, by = 'decile')
# reorder columns: decile, method, 75%, n_family, n_ref_branches, n_tot
summary_table = summary_table[,c("decile", "method", "75%", "n_families", "n_ref_branches", "n_tot")]
# rename columns 
colnames(summary_table) = c('Decile', 'Method', '#cols for 75% ref branches', '#families', '#reference branches', 'Number of columns')
# write table
tables = paste(source_data, '/../tables', sep = '')
write.table(summary_table, file=paste(tables, paste("Titration_deciles.csv", sep=''), sep = "/"), append=F,quote=F, sep = ",", eol = "\n", na = "NA", dec = ".",row.names =F, col.names =T)




