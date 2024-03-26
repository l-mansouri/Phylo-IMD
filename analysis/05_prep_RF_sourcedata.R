library(phangorn)


# ! PIPELINE OUTPUT FOLDER, CHANGE ACCORDINGLY
output_folder = '/home/luisasantus/Desktop/crg_cluster/NF_draft/'
# ! SOURCE DATA FOLDER, CHANGE ACCORDINGLY
source_data = '/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data/'

# -----------------------------------------------------------------------------
# Overview of the script: 
# read in trees and calculate RF distance
# -----------------------------------------------------------------------------


fl=read.table(paste("source_data",'list_of_families_with_all_rep_in_3d', sep = "/"))[,1]
al="mTMalign"
tr="untrimmed"
RF3d=c()
RFtm=c()
setwd(output_folder)
for (fam in fl){
  #importing ML tree
  trml = read.tree(paste(al, '_ML_', tr,'_trees/',fam, '_', al, '_ML_', tr, '.nwk',sep=''))
  tr3d = read.tree(paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_ME_',tr,'.nwk', sep=''))
  #importing tm score matrix and tree if it's computed
  trtm = read.tree(paste('mTMalign_TM_ME_trees/',fam,'_mTMalign_TM_ME_untrimmed.nwk',sep=''))
  rf3dml=RF.dist(tr3d, trml, normalize=T)
  rftmml=RF.dist(trtm, trml, normalize=T)
  RF3d=c(RF3d,rf3dml)
  RFtm=c(RFtm, rftmml)
}
df=data.frame(RF3d, RFtm)

write.table(df, file = paste(source_data, "RF_table.csv", sep = "/"), sep = ",", qmethod = "double", row.names = FALSE)


