library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)

# ---------------------------- MAIN ------------------------------------------------------------------
setwd('/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/')
mtm=read.table('source_data/mTMalign_ME_untrimmed_saturation_correlations.txt', header=T)

# mtmALIGN UNTRIMMED
al="mTMalign"
tr="untrimmed"

# Square everything 
mtm = mtm^2

median_seq = round(median(mtm$COR1dme),3)
median_imd = round(median(mtm$COR3dme),3)
median_tm  = round(median(mtm$CORtmme),3)
label_seq  = paste("ME \n (median: ",median_seq, ")", sep = "")
label_imd  = paste("IMD \n (median: ",median_imd, ")", sep = "")
label_tm  = paste("TM \n (median: ",median_tm, ")", sep = "")
label_seq = "ME"
label_tm = "TM"
label_imd = "IMD"
df111=data.frame(correlations=c(mtm$COR1dme, mtm$CORtmme, mtm$COR3dme),
    type=factor(rep(c(label_seq, label_tm, label_imd), each=length(mtm[,1])), levels=c(label_seq, label_tm, label_imd)))

palette = c("#27A592", "#2191FB", 158 )
palette = c("#A331D8","#27A592", "#2191FB" )
#palette = c(158,"#27A592", "#2191FB" )

col_seq = palette[1]
col_tm  = palette[2]
col_imd = palette[3]

p_sat1=ggplot(df111, aes(x=type, y=correlations, fill=type, col = type, alpha = 0.2))+
    stat_boxplot(geom = "errorbar",
               width = 0.15, linetype = "solid") + 
    geom_boxplot(outlier.size=0.1, lwd=0.6, linetype = "solid" )+
    xlab('')+
    ylab('R² (input distances, patristic distances)')+
    annotate('text', x=df111$type[1], y=1.1, label=label_seq, vjust=1.7, col = col_seq, size = 6)+
    annotate('text', x=df111$type[513], y=1.1, label=label_tm, vjust=1.7, col = col_tm, size = 6)+
    annotate('text', x=df111$type[1025], y=1.1, label=label_imd, vjust=1.7, col = col_imd, size = 6)+
    ylim(0,1.19)+
    scale_fill_manual(values=palette)+
    scale_color_manual(values=palette)+
    theme_minimal()+
    theme(legend.position = "none",
          axis.text = element_text(color = "black", size = 10),
          axis.text.x = element_blank(),    # Hide axis text
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 13),
          axis.ticks.x = element_blank(),    # Hide axis ticks
          panel.grid = element_blank(), 
          plot.background = element_rect(fill = "white"),
          axis.ticks.y = element_line(color = "black")) +
    #geom_errorbar(linetype = "dashed")+
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1), linetype = "solid", color = "black")+
    annotate('text', x=df111$type[1], y=0.4, label=( paste("median =", median_seq)), vjust=1.7, col = "black", size = 3.5)+
    annotate('text', x=df111$type[513], y=0.53, label=(paste("median =", median_tm)), vjust=1.7, col = "black", size = 3.5)+
    annotate('text', x=df111$type[1025], y=0.53, label=(paste("median =", median_imd)), vjust=1.7, col = "black", size = 3.5)+
  theme(panel.background = element_rect(fill = "transparent", color = NA))+ theme(plot.background = element_rect(color = NA))
ggsave('plots/main/Fig1_correlation_input_patristic_mtmuntrimmed_508.png', plot=p_sat1, width = 6, height = 7, dpi = 300)

mtm

# mtm df columns COR1dme, CORtmme, COR3dme and rename as Seq-ME, TM-ME, IMD-ME
mtm = mtm[,c("COR1dme", "CORtmme", "COR3dme")]
colnames(mtm) = c("SeqME", "TMME", "IMDME")

# Paired wilcoxon test - test for significance, wilcoxon paired test for IMD higher than Seq
wilcox.test(mtm$IMDME, mtm$SeqME, paired = TRUE, alternative = 'greater')
# Paired wilcoxon test - test for significance, wilcoxon paired test for TM higher than Seq
wilcox.test(mtm$TMME, mtm$SeqME, paired = TRUE, alternative = 'greater')

# Compute median of each column
apply(mtm, 2, median)

# ---------------------------- SUPPL ------------------------------------------------------------------


plot_correlation_patristic <- function(al,tr, tag = ""){
  mtm=read.table(paste('source_data/',al,'_',tr,'_saturation_correlations_self_ML.txt', sep = ""), header = T)

  if( al == "tcoffee" ){
    al = "T-Coffee"
  }else if( al == "sap_tmalign"){
    al = "3D-Coffee"
  }
  
  title = paste(al,tr, sep  = " ")
  
  median_seq = round(median(mtm$COR1dme),3)
  median_imd = round(median(mtm$COR3dme),3)
  label_seq  = paste("ME \n (median: ",median_seq, ")", sep = "")
  label_imd  = paste("IMD \n (median: ",median_imd, ")", sep = "")
  
  df111=data.frame(correlations=c(mtm$COR1dme, mtm$COR3dme),
                   type=factor(rep(c(label_seq, label_imd), each=length(mtm[,1])), levels=c(label_seq, label_imd)))
  
  
  p_sat1=ggplot(df111, aes(x=type, y=correlations, fill=type, col = type, alpha = 0.2))+
    stat_boxplot(geom = "errorbar",
                 width = 0.15, linetype = "solid") + 
    geom_boxplot(outlier.size=0.1, lwd=0.6, linetype = "solid" )+
    xlab('')+
    ylab('R² (input distances, patristic distances)')+
    annotate('text', x=df111$type[1], y=1.1, label=label_seq, vjust=1, col = col_seq, size = 2.8)+
    annotate('text', x=df111$type[513], y=1.1, label=label_imd, vjust=1, col = col_imd, size = 2.8)+
    ylim(0,1.19)+
    scale_fill_manual(values=palette[c(1,3)])+
    scale_color_manual(values=palette[c(1,3)])+
    theme_minimal()+
    theme(legend.position = "none",
          axis.text = element_text(color = "black", size = 8),
          axis.text.x = element_blank(),    # Hide axis text
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.ticks.x = element_blank(),    # Hide axis ticks
          panel.grid = element_blank(), 
          plot.background = element_rect(fill = "white"),
          axis.ticks.y = element_line(color = "black")) +
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1), linetype = "solid", color = "black")+
    labs(title = title, tag = bquote(bold(.(tag))))+
    theme(plot.title = element_text(hjust = 0.6, vjust = -12, size = 10))+
    theme(panel.background = element_rect(fill = "transparent", color = NA))
  
  
  return(p_sat1)
}

p1 <- plot_correlation_patristic('mTMalign',"trimmed", "A") 
p2 <- plot_correlation_patristic('sap_tmalign',"untrimmed", "B")
p3 <- plot_correlation_patristic('sap_tmalign',"trimmed", "C")
p4 <- plot_correlation_patristic('tcoffee',"untrimmed", "D")
p5 <- plot_correlation_patristic('tcoffee',"trimmed", "E")


panel <- grid.arrange(p1 + theme(plot.background = element_rect(fill = "transparent", color = NA)),
                               p2 + theme(plot.background = element_rect(fill = "transparent", color = NA)),
                               p3 + theme(plot.background = element_rect(fill = "transparent", color = NA)),
                               p4 + theme(plot.background = element_rect(fill = "transparent", color = NA)),
                               p5 + theme(plot.background = element_rect(fill = "transparent", color = NA)),
                               ncol = 3)

ggsave(paste('plots/suppl/S2_correlation_input_patristic.png', sep = ""), plot=panel, width = 9, height = 10, dpi = 300)



