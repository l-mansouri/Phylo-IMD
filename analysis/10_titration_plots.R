library(ggplot2)
library(patchwork)
# ------------------------------  MAIN -------------------------------------------------------

theme_set(theme_bw())
theme_set(theme_bw(base_size = 18))
theme_update(legend.title=element_blank())

# ----------------------------------------
#           FIGURE 4A
# ----------------------------------------

plot_titration <- function(df, labels, palette =  c("#f8766d", "#00ba39", "#619cff")){
  
  colnames(df)=c('ncol','avg1d', 'sd1d', 'avg3d', 'sd3d', 'avgML', 'sdML', 'x')
  dfml=data.frame(n_o_col=c(df$ncol, df$ncol, df$ncol),
                  avg_RF=c(df$avg1d, df$avg3d, df$avgML),
                  sd_RF=c(df$sd1d, df$sd3d, df$sdML),
                  type=factor(rep(labels, each=as.numeric(length(df[,1]))), levels=labels)
  )
  my_palette <- palette

  # Plot using custom palette
  pexp=ggplot(dfml, aes(x=dfml[,1], y=dfml[,2], color=type), alpha=0.2)+
    geom_point()+
    geom_line()+
    xlab('Number of columns included')+
    ylab('Fraction of reference branches present')+
    geom_ribbon(aes(ymin = dfml[,2] - dfml[,3], ymax = dfml[,2] + dfml[,3], fill = type), color = NA, alpha = 0.2)+
    scale_color_manual(values = my_palette)+
    scale_fill_manual(values = my_palette)+
    xlim(c(0,201))+
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), limits = c(0,1.03))
  pexp=pexp+geom_hline(yintercept=0.75, linetype="dashed", color = "grey", linewidth=1)
  return(pexp)  
}


# Fig. 4A ( all families, experimental structures )
df_exp_full=read.table ('source_data/titration_every_5_reference_branches_ME+3d+ML.txt')
pexp= plot_titration(df_exp_full,c("ME", "IMD", "ML")) 
#ggsave('plots/main/Fig4A_fraction_of_ref_branches_ME+3d+ML_avg_replicates.png', plot=pexp,  width = 7, height = 5.7, dpi = 300, units = "in")

# Fig. 4B (families with AF2 )
df_af2=read.table ('source_data/af2_pred_pdb_titration_every_5_reference_branches_ME+3d+ML.txt')
#paf2_only= plot_titration(df_af2,c("AF2 ME", "AF2 IMD", "AF2 ML"), palette = c("#910f06","#0e8a00", "#040052")) 
paf2_only= plot_titration(df_af2,c("AF2 ME", "AF2 IMD", "AF2 ML")) 

#ggsave('plots/main/Fig4B_AF2only_comparison_af2-exp_pdb_fraction_of_ref_branches_ME+3d+ML_avg_replicates.png', plot=paf2_only, width = 7, height = 5.7, dpi = 300, units = "in")

# MAIN FIGURE 4 PANEL 
panel <- (pexp+paf2_only)+ plot_annotation(tag_levels = 'A')+ plot_layout(ncol = 2)
ggsave(paste('plots/',"main",'/Fig4_titration.png', sep = ""), plot=panel, width = 15, height =6, dpi = 300)


# SUPPLEMENTARY 
df_3dcoffe = read.table('source_data/titration_every_5_reference_branches_ME+3d+ML_3DCOFFEE.txt')
p3dcoffee= plot_titration(df_3dcoffe,c("ME", "IMD", "ML")) 
p3dcoffee = p3dcoffee+ggtitle("Titration on 3D-Coffee MSAs")
#ggsave('plots/main/titration_3dcoffee.png', plot=p3dcoffee, width = 7, height = 5.7, dpi = 300, units = "in")



# PLOT BOTH (exp and af2)
df1=read.table ('source_data/exp_pdb_titration_every_5_reference_branches_ME+3d+ML.txt')
colnames(df1)=c('ncol','avg1d', 'sd1d', 'avg3d', 'sd3d', 'avgML', 'sdML', 'x')
colnames(df_af2)=c('ncol','avg1d', 'sd1d', 'avg3d', 'sd3d', 'avgML', 'sdML', 'x')

DF=data.frame(n_o_col=c(df1$ncol, df1$ncol, df1$ncol, df_af2$ncol, df_af2$ncol, df_af2$ncol),
              avg_RF=c(df1$avg1d, df1$avg3d, df1$avgML, df_af2$avg1d, df_af2$avg3d, df_af2$avgML),
              sd_RF=c(df1$sd1d, df1$sd3d, df1$sdML, df_af2$sd1d, df_af2$sd3d, df_af2$sdML),
              type=factor(rep(c("exp ME", "exp IMD", "exp ML", "AF2 ME", "AF2 IMD", "AF2 ML"), each=as.numeric(length(df1[,1]))), levels=c("exp ME",  "AF2 ME", "exp IMD","AF2 IMD", "exp ML", "AF2 ML"))
)
my_palette <- c("#f8766d", "#910f06", "#00ba39", "#0e8a00", "#619cff", "#040052")

paf2=ggplot(DF, aes(x=DF[,1], y=DF[,2], color=type), alpha=0.2)+geom_point()+geom_line()+xlab('Number of columns included')+ylab('Fraction of reference branches present')+scale_color_manual(values = my_palette)
paf2=paf2+geom_hline(yintercept=0.75, linetype="dashed", color = "grey", size=1)
paf2 = paf2+geom_ribbon(aes(ymin = DF[,2] - DF[,3], ymax = DF[,2] + DF[,3], fill = type), color = NA, alpha = 0.2)+
  scale_color_manual(values = my_palette)+
  scale_fill_manual(values = my_palette)+
  xlim(c(0,201))+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), limits = c(0,1.03))

# make the figure bigger
# ggsave('plots/suppl/Titration_comparison_af2-exp_pdb_fraction_of_ref_branches_ME+3d+ML_avg_replicates.png', plot=paf2, width = 9.5, height = 8, dpi = 300, units = "in")


#--------------------------------
#   panel for supplementary 
#--------------------------------
panel_suppl = (p3dcoffee+paf2)+ plot_annotation(tag_levels = 'A')+ plot_layout(ncol = 2)
ggsave(paste('plots/',"suppl",'/S7_titration.png', sep = ""), plot=panel_suppl, width = 15, height =6, dpi = 300)






