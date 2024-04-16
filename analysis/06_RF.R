library(ggplot2)
library(cowplot)
library(RColorBrewer)

# -----------------------------------------------------------------------------
# Overview of the script: 
# read in trees and calculate RF distance
# -----------------------------------------------------------------------------

df = read.table("source_data/RF_table_rooted.csv", header = TRUE, sep = ",")
okRFme=df[df$RF3d<df$RFtm,1]
sameRFme=df[df$RF3d==df$RFtm,1]
df$delta=(df$RF3d - df$RFtm)
p1=ggplot(df, aes(x=RF3d, y=RFtm))+geom_abline(color='lightgrey')+geom_point(size = 1)+xlab('RF(IMD-ME, Seq-ML)')+ylab('RF(TM-ME, Seq-ML)')+
  theme_bw()+annotate("label", label=length(okRFme), x=0.05, y=1, color='springgreen4')+
  annotate("label", label=508-length(sameRFme)-length(okRFme), y=0.01, x=1, color='red')+
  annotate("label", label=length(sameRFme), x=1.05, y=1.05, color='gray40')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 10))


source('SIGN.test.R')
SIGN.test(df$delta, md=0, alternative='less')

p1=p1+annotate('text', x=0.85, y=0.2, label='SIGN test\np-value = 0.008')
ggsave('plots/main/Fig2_mTMalign_RF_3d-TM_ML_signtest.png', p1, width = 5, height = 5) 
