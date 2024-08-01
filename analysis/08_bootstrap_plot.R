library(grid)
library(gridExtra)
library(RColorBrewer)
library(patchwork)

fl=read.table('source_data/list_of_families_with_all_rep_in_3d')[,1]


# ---------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------

scatter_plot <- function(df, title = "", tag = list("A", "B")){
    df = read.table(df,   header = TRUE, sep = ",")

    C3dme=round(cor(df$average_1d, df$average_3d),2)
    C3dml=round(cor(df$average_ML, df$average_3d),2)
    CrfME=round(cor(df$average_3d,df$RF_3d_ME ),2)
    CrfML=round(cor(df$average_3d,df$RF_3d_ML ),2)

    cRFme=c()
    cRFml=c()
    xpos=c()
    for (th in 1:100){
        crfme=c()
        crfml=c()
        for (p in 1:length(df[,1])){
            r =as.numeric(df[p,])
            if (r[2]>=th && r[2]<th+10){
                crfme=c(crfme,r[5])
                crfml=c(crfml, r[6])
            }
        }
        xpos=c(xpos,th)
        cRFme=c(cRFme, mean(crfme))
        cRFml=c(cRFml, mean(crfml))
    }

    df_rf=data.frame(xpos, cRFme,cRFml)
    df_rf1=data.frame(xpos, cRFme)
    df_rf2=data.frame(xpos, cRFml)

    coeff=100

    p1<- ggplot(df, aes(y=average_1d, x=average_3d, color=RF_3d_ME) ) + geom_point() + scale_color_gradient(low="red", high="blue")+
        ylab('ME average bootstrap support value')+
        xlab('IMD average bootstrap support value')+
        labs(color='RF') +geom_abline(color='gray')+
        annotate('text', x=35, y=95, label=paste('R = ', C3dme,sep='')) +
        xlim(c(0,100))+ylim(c(0,100))+theme_light()+
        theme(axis.text = element_text( size = 12, color = "black"))+
        theme(axis.title = element_text( size = 14, color = "black"))
            

        
    df_rf1=na.omit(df_rf1)
    
    p12=p1+ geom_line(data=df_rf1, aes(x=xpos, y=cRFme*coeff), color='black', size=1.2, alpha=0.7) + 
        scale_y_continuous(
        # Features of the first axis
        #name = paste('average bootstrap ',al,'_Seq-ME', sep=''), 
        # Add a second axis and specify its features
        sec.axis = sec_axis(~./coeff, name="Average RF"),
        limits=c(0,100)
    ) + labs(title = title, tag = bquote(bold(.(tag[[1]]))))+
        theme(axis.title = element_text( size = 14))+ theme(plot.title = element_text(hjust = 0.5))




    p2<- ggplot(df, aes(y=average_ML, x=average_3d, color=RF_3d_ML) ) + geom_point() + scale_color_gradient(low="red", high="blue")+
        xlab('IMD average bootstrap support value')+
        ylab('ML average bootstrap support value')+
        labs(color='RF') +geom_abline(color='gray')+
        annotate('text', x=33, y=95, label=paste('R = ', C3dml,sep=''))+
        xlim(c(0,100))+ylim(c(0,100))+theme_light()+
        theme(axis.text = element_text( size = 12, color = "black"))+
        theme(axis.title = element_text( size = 14, color = "black"))+
        labs(title = title, tag = bquote(bold(.(tag[[2]]))))
        theme(axis.title = element_text( size = 14))+ theme(plot.title = element_text(hjust = 0.5))



    df_rf2=na.omit(df_rf2)

    p23=p2+ geom_line(data=df_rf2, aes(x=xpos, y=cRFml*coeff), color='black', size=1.2, alpha=0.7) + 
    scale_y_continuous(
    # Features of the first axis
    #name = paste('average bootstrap ',al,'_Seq-ML', sep=''), 
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeff, name="Average RF"),
    limits=c(0,100)
    )

    return(list(p12, p23))
}


# --------------------- Second part of pane: lineplot ---------------------------------------------------------------------------


plot_lines <- function(al,tr, folder, tag = list("","")){
  DF= read.table(paste('source_data/',al,'_',tr,'counts_of_found_in_ML_and_ME_with_moving_IMD_BS_thr.txt', sep=''), header = T)
  DF["TTB"]  = max(DF["TB"])
  OKB_ML = DF[,"OKB_ML"]
  OKB_ME = DF[,"OKB_ME"]
  TTB    = DF[,"TTB"]
  TB = DF[,"TB"]
  
  if( al == "tcoffee" ){
    al = "T-Coffee"
  }else if( al == "sap_tmalign"){
    al = "3D-Coffee"
    print(al)
  }
  
  title = ""
  if(folder == "suppl" ){
    title = paste(al,tr,sep = " ")
  }
  
  
  dfML=data.frame(pos=c(1:100, 1:100),
                  frac_ok=c((OKB_ML/TB)*100, (TB/TTB)*100),
                  type=factor(rep(c('found in ML', 'fraction over total\nnumber of\nbranches'), each=100), levels=c('found in ML', 'fraction over total\nnumber of\nbranches'))
  )
  dfME=data.frame(pos=c(1:100, 1:100),
                  frac_ok=c((OKB_ME/TB)*100, (TB/TTB)*100),
                  type=factor(rep(c('found in ME', 'fraction over total\nnumber of\nbranches'), each=100), levels=c('found in ME', 'fraction over total\nnumber of\nbranches'))
  )
  color1 = "salmon"
  color2 = "#00BFC4"
  palette <-c(color1, color2)
  
  PML=ggplot(dfML, aes(x=pos, y=frac_ok, color=type))+geom_line(size=1)+xlab('IMD average bootstrap suport value')+labs(title = title)+scale_y_continuous(
    name = 'Branches occurring in ML (%)', 
    sec.axis = sec_axis(~., name="Considered branches (%)"),
    limits=c(0,100),
    breaks = c(0,25,50,75,100)
  )+theme_light()+
    scale_color_manual(values = palette)+
    theme(axis.title.y.right = element_text(color = color2, angle = 90, face = "bold"),
          axis.title.y.left = element_text(color = color1, face = "bold"),
          axis.ticks.x = element_blank(),
          legend.position = "none")+
    scale_x_continuous(breaks = c(0,25,50,75,100))+
    theme(axis.text = element_text( size = 12, color = "black"))+
    theme(axis.title = element_text( size = 14))+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text = element_text(color = "black", size = 10))+
    labs(title = title, tag = bquote(bold(.(tag[[2]]))))+
    theme(plot.margin = margin(1, 1, 1, 1, "cm"))
    
  
  PME=ggplot(dfME, aes(x=pos, y=frac_ok, color=type))+geom_line(size=1)+xlab('IMD average bootstrap support value')+labs(title = title)+scale_y_continuous(
    name = 'Branches occurring in ME (%)', 
    sec.axis = sec_axis(~., name="Considered branches (%)"),
    limits=c(0,100),
    breaks = c(0,25,50,75,100)
  )+theme_light()+
    theme(axis.text = element_text(color = "black", size = 10))+
    scale_color_manual(values = palette)+
    theme(axis.title.y.right = element_text(color = color2, angle = 90, face = "bold"),
          axis.title.y.left = element_text(color = color1, face = "bold"),
          axis.ticks.x = element_blank(),
          legend.position = "none")+
    scale_x_continuous(breaks = c(0,25,50,75,100))+
    theme(axis.text = element_text( size = 12, color = "black"))+
    theme(axis.title = element_text( size = 14))+ theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text = element_text(color = "black", size = 10))+
    theme(panel.background = element_rect(fill = "transparent", color = NA))+
    labs(title = title, tag = bquote(bold(.(tag[[1]]))))+
    theme(plot.margin = margin(1, 1, 1, 1, "cm"))
  
  
  return(list(PME,PML))
}

# ---------------------------------
# MAIN
# --------------------------------- 

main <- scatter_plot("source_data/mTMalign_untrimmed_bootstrap_table.csv",tag = list("A", "B"))
pA <- main[[1]]
pB <- main[[2]]


main <- plot_lines("mTMalign", "untrimmed", "main", tag = list("C", "D"))
pC <- main[[1]]
pD <- main[[2]]

# Arrange the plots with spacing
panel <- (pA+pB)/(pC+pD)

ggsave(paste('plots/',"main",'/Fig3_bootstrap.png', sep = ""), plot=panel, width = 14, height =12, dpi = 300)


# ---------------------------------
# SUPPLEMENTARY
# ---------------------------------

# --------------- SCATTERPLOTS ---------------------

# UNTRIMMED
pAB <- scatter_plot("source_data/sap_tmalign_untrimmed_bootstrap_table.csv", title = "3D-Coffee untrimmed", tag = list("A", "B"))
pa <- pAB[[1]]
pb <- pAB[[2]]
pCD <- scatter_plot("source_data/tcoffee_untrimmed_bootstrap_table.csv", title = "T-Coffee untrimmed", tag = list("C", "D"))
pc <- pCD[[1]]
pd <- pCD[[2]]

panel <- (pa+pb)/(pc+pd)
ggsave(paste('plots/',"suppl",'/S8_scatter_bootstrap_untrimmed.png', sep = ""), plot=panel, width = 14 , height =12)


# TRIMMED
pAB <- scatter_plot("source_data/mTMalign_trimmed_bootstrap_table.csv", title = "mTMalign trimmed", tag = list("A", "B"))
pa <- pAB[[1]]
pb <- pAB[[2]]
pCD <- scatter_plot("source_data/sap_tmalign_trimmed_bootstrap_table.csv", title = "3D-Coffee trimmed", tag = list("C", "D"))
pc <- pCD[[1]]
pd <- pCD[[2]]
pEF <- scatter_plot("source_data/tcoffee_trimmed_bootstrap_table.csv", title = "T-Coffee trimmed", tag = list("E", "F"))
pe <- pEF[[1]]
pf <- pEF[[2]]

panel <- (pa+pb)/(pc+pd)/(pe+pf)
ggsave(paste('plots/',"suppl",'/S9_scatter_bootstrap_trimmed.png', sep = ""), plot=panel, width = 15 , height =18)


#---------------- LINEPLOTS ---------------------
# UNTRIMMED
pAB <- plot_lines("sap_tmalign", "untrimmed", "suppl", list("A", "B"))
pa <- pAB[[1]]
pb <- pAB[[2]]
pCD <- plot_lines("tcoffee", "untrimmed", "suppl", list("C", "D"))
pc <- pCD[[1]]
pd <- pCD[[2]]

panel <- (pa+pb)/(pc+pd)
ggsave(paste('plots/',"suppl",'/S5_lineplot_bootstrap_untrimmed.png', sep = ""), plot=panel, width = 9.5 , height =9)

# TRIMMED
pAB_t <- plot_lines("mTMalign", "trimmed", "suppl", list("A", "B"))
pa_t <- pAB_t[[1]]
pb_t <- pAB_t[[2]]
pCD_t <- plot_lines("sap_tmalign", "trimmed", "suppl", list("C", "D"))
pc_t <- pCD_t[[1]]
pd_t <- pCD_t[[2]]
pEF_t <- plot_lines("tcoffee", "trimmed", "suppl", list("E", "F"))
pe_t <- pEF_t[[1]]
pf_t <- pEF_t[[2]]

panel <- (pa_t+pb_t)/(pc_t+pd_t)/(pe_t+pf_t)
ggsave(paste('plots/',"suppl",'/S6_lineplot_bootstrap_trimmed.png', sep = ""), plot=panel, width = 9.5 , height =13, dpi = 300)


