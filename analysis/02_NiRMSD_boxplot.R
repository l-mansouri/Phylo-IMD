library(ggplot2)

DF2 = read.table("source_data/S1_NiRMSD_table.csv", header = TRUE, sep = ",")
DF2$NiRMSD <- as.numeric(DF2$NiRMSD)
DF2$type <- factor(DF2$type, levels = c('T-Coffee', '3D-Coffee', 'mTM-align'))

df_tc = DF2[DF2$type == "T-Coffee",]
df_st = DF2[DF2$type == "3D-Coffee",]
df_mt = DF2[DF2$type == "mTM-align",]


palette = c("#4421af", "darkgrey","#ef9b20")

p1=ggplot(DF2, aes(x=type, y=NiRMSD, fill=type, col = type, alpha = 0.6))+
  stat_boxplot(geom = "errorbar",
               width = 0.15, linetype = "solid") + 
  geom_boxplot(outlier.size=0.5, lwd=0.6, linetype = "solid", notch = T )+
  xlab('')+
  ylab('Average NiRMSD per family')+
  theme_light()+
  theme(legend.position="None")+
  annotate('text', x=DF2$type[1], y=median(df_tc[,1]), label=round(median(df_tc[,1]),2), vjust=1.2, size = 5)+
  annotate('text', x=DF2$type[length(df_tc[,1])+1], y=median(df_st[,1]), label=round(median(df_st[,1]),2), vjust=1.2, size= 5)+
  annotate('text', x=DF2$type[(length(df_tc[,1])*2)+1], y=median(df_mt[,1]), label=round(median(df_mt[,1]),2), vjust=1.2, size = 5)+
  theme(axis.text = element_text( color = "black", size = 13), axis.title = element_text( size = 15 ))+
  scale_fill_manual(values=palette)+
  scale_color_manual(values=palette)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black", size = 13),
        axis.text.x = element_text(size = 15),    # Hide axis text
        axis.title.y = element_text(size = 15),
        axis.ticks.x = element_blank(),    # Hide axis ticks
        panel.grid = element_blank(), 
        plot.background = element_rect(fill = "white"),
        axis.ticks.y = element_line(color = "black")) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 2.5), linetype = "solid", color = "black")+
  theme(panel.background = element_rect(fill = "transparent", color = NA))+ theme(plot.background = element_rect(color = NA))

ggsave('plots/suppl/S1_NiRMSD_boxplot_508.png', plot=p1, width = 7, height = 7)
