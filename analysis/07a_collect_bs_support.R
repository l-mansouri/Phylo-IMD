

# ! PIPELINE OUTPUT FOLDER, CHANGE ACCORDINGLY
output_folder = '/home/luisasantus/Desktop/crg_cluster/NF_draft/'
setwd(output_folder)

# -----------------------------------------------------------------------------
# Overview of the script: 
# read in replicates
# calculate bootstrap support for each node
# write bootstrap support to tree
# -----------------------------------------------------------------------------


library(phangorn)
aligners=c('mTMalign', 'sap_tmalign', 'tcoffee')
trimming=c('untrimmed', 'trimmed')

for (al in aligners){
    if (al=='mTMalign'){
        fl=read.table('mTMalign_list')[,1]
    }else{
        fl=read.table('list_of_families')[,1]
    }
    for (tr in trimming){
        if (al=='mTMalign' && tr=='trimmed'){
            fl=fl[-102]
            fl=fl[-222]
            fl=fl[-491]
        }
        if (al=='tcoffee' && tr=='trimmed'){
            fl=fl[-4]
        }
        for (fam in fl){
            #importing 3d tree
            if (al=='mTMalign'){
                tree3d = read.tree(paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_ME_',tr,'.nwk', sep=''))
            } else{
                tree3d = read.tree(paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_',tr,'.nwk', sep=''))
            }
            bootrees3d=read.tree(paste(al,'_3d_ME_',tr,'_trees/replicates/', fam,'_bootrees', sep=''))
            b3d=prop.clades(tree3d, bootrees3d)
            b3d[is.na(b3d)]<-0
            tree3d$node.label=b3d
            tree3d$node.label[1]=''
            if (al=='mTMalign'){
                #write.tree(tree3d, file=paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_ME_',tr,'.nwk', sep=''))
            } else{
                #write.tree(tree3d, file=paste(al,'_3d_ME_',tr,'_trees/',fam,'_',al,'_3d_',tr,'.nwk', sep=''))
            }
        }
    }
}