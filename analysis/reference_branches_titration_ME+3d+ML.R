#!/usr/bin/env Rscript
library('phangorn')
library('geiger')
args = commandArgs(trailingOnly=TRUE)


fam=args[1]


no_col=c()
total_shared=c()
agreement1d1d=c()
agreement3d3d=c()
agreementmlml=c()
total_total=c()

for (fr in seq(from=10, to=200, by=5)){
  total_tot=0
  total_sh=0
  ag1d1d=0
  ag3d3d=0
  agmlml=0
  no_col=c(no_col,fr)
  for (rep in 0:9){
    #reading the trees computed on the 200 columns
    tree1d=read.tree(paste('exp_pdb_results_SUBSET49/ME_trees/',fam,'_mTMalign_untrimmed_random_msa_replicate.',rep,'_with_200_columns.nwk', sep=''))
    tree3d=read.tree(paste('exp_pdb_results_SUBSET49/IMD_trees/',fam,'_mTMalign_untrimmed.random_column_pairs_replicate.',rep,'_with_200_columns.txt.matrices_fastme_tree.nwk', sep=''))
    treeml=read.tree(paste('exp_pdb_results_SUBSET49/ML_trees/',fam,'_mTMalign_untrimmed_random_msa_replicate.',rep,'_with_200_columns.ph.treefile', sep=''))
    #computing the shared branches: using prop.clades with 2 trees the possible values are 0,1,2 and we only need the 2; that have to be transformed into 1s for the sake of this script hence the modification of the vector
    v1d=prop.clades(tree1d, tree3d, treeml)
    v1d[is.na(v1d)] <- 0
    v1d[v1d==1] <- 0
    v1d[v1d==2] <- 1
    v1d[1]=0
    v3d=prop.clades(tree3d, tree1d, treeml)
    v3d[is.na(v3d)] <- 0
    v3d[v3d==1] <- 0
    v3d[v3d==2] <- 1
    v3d[1]=0
    vml=prop.clades(treeml, tree1d, tree3d)
    vml[is.na(vml)] <- 0
    vml[vml==1] <- 0
    vml[vml==2] <- 1
    vml[1]=0
    total_sh=total_sh+sum(v1d) #number of shared branches
    total_tot=total_tot+(length(v1d)-1)
    #importing the trees on the fraction
    ftree1d=read.tree(paste('exp_pdb_results_SUBSET49/ME_trees/',fam,'_mTMalign_untrimmed_random_msa_replicate.',rep,'_with_',fr,'_columns.nwk', sep=''))
    ftree3d=read.tree(paste('exp_pdb_results_SUBSET49/IMD_trees/',fam,'_mTMalign_untrimmed.random_column_pairs_replicate.',rep,'_with_',fr,'_columns.txt.matrices_fastme_tree.nwk', sep=''))
    ftreeml=read.tree(paste('exp_pdb_results_SUBSET49/ML_trees/',fam,'_mTMalign_untrimmed_random_msa_replicate.',rep,'_with_',fr,'_columns.ph.treefile', sep=''))
    #computing the shared br between the fraction tree and the full(200 col) tree
    pc1d1d=prop.clades(tree1d, ftree1d)
    pc1d1d[is.na(pc1d1d)] <- 0
    pc1d1d[1]=0
    pc3d3d=prop.clades(tree3d, ftree3d)
    pc3d3d[is.na(pc3d3d)] <- 0
    pc3d3d[1]=0
    pcmlml=prop.clades(treeml, ftreeml)
    pcmlml[is.na(pcmlml)] <- 0
    pcmlml[1]=0
    #computing the agreement
    ag1d1d=ag1d1d+sum(pc1d1d*v1d)
    ag3d3d=ag3d3d+sum(pc3d3d*v3d)
    agmlml=agmlml+sum(pcmlml*vml)
  }
  #computing the fraction of shared branches
  ag1d1d=ag1d1d/total_sh
  ag3d3d=ag3d3d/total_sh
  agmlml=agmlml/total_sh
  #collecting the results
  agreement1d1d=c(agreement1d1d, ag1d1d)
  agreement3d3d=c(agreement3d3d, ag3d3d)
  agreementmlml=c(agreementmlml, agmlml)
  total_shared=c(total_shared, total_sh)
  total_total=c(total_total, total_tot)
}


df=data.frame(no_col, agreement1d1d, agreement3d3d , agreementmlml, total_shared, total_total)
#df[is.na(df[,])]=0
# create directory if it doesn't exist
dir.create('exp_pdb_results_SUBSET49/ref_br_ME+3d+ML', showWarnings = FALSE)
write.table(df, file=paste('exp_pdb_results_SUBSET49/ref_br_ME+3d+ML/', fam,'_avg_fr_ref_br_per_fraction_shared_1d_3d_ML', sep=''), append=F,quote=F, sep = " ", eol = "\n", na = "0", dec = ".",row.names = F, col.names =F)
