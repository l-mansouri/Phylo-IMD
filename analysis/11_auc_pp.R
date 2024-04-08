# --------------------------------------------------
#  Compute the references 
# --------------------------------------------------
#!/usr/bin/env Rscript
library('phangorn')
library('geiger')
source("utils_auc.R")


source_data = "/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data"
# OUTPUT DIRECTORY OF SPLIT FILES 
output_dir = "/home/luisasantus/Desktop/crg_cluster/newphylo/NF_TMalign/mTMalign/titration_every_5_bootstrap_200_columns/split_files"


# FAMILIES 
families =read.table(paste(source_data,'list_of_families_for_titration', sep = "/"))[,1]

# BOOTSTRAP THRESHOLD
bs_threshold = 80 



#--------------------------------------------------
# --------------------------------------------------
#                     Main
# --------------------------------------------------
# --------------------------------------------------

# --------------------------------------------------
#  SUPPLEMENTARY TABLE 3 
# --------------------------------------------------


compute_references <- function(families, dir, dir_2, dir_3, method, bs_threshold){
    ref_count        = 0
    count_trees_used = 0 
    tot_branches     = 0
    depths           = c() 
    relative_depths  = c()
    # FOR EACH FAMILY 
    for (fam in families ){
        #-------------------------
        # 1. Get the splits
        #-------------------------
        # If only one split file is provided
        split_file = parse_splits(fam,dir, output_dir)
        # If two split files are provided
        if (dir_2 != ""){
            split_file_2 = parse_splits(fam,dir_2, output_dir)
            splits_combined_bs = combine_bs(split_file, split_file_2, method)
            # If three split files are provided
            if(dir_3 != ""){
                split_file_3 = parse_splits(fam,dir_3, output_dir)
                splits_combined_bs = combine_bs(splits_combined_bs, split_file_3, method)
            }   
        }
        else {
            splits_combined_bs = split_file
        }
        references  = get_splits(splits_combined_bs, bs_threshold)
        #-------------------------
        # 2. Get the infos
        #-------------------------
        if(length(references) > 0){
            ref_count = ref_count + length(references)
            total_branches_family = length(get_splits(splits_combined_bs, 0))
            nseq = total_branches_family + 3
            tot_branches = tot_branches + total_branches_family
            fam_depths = get_depths(splits_combined_bs, bs_threshold)
            depths = c(depths, fam_depths)
            relative_depths = c(relative_depths, ((fam_depths/(nseq))*2)) 
            count_trees_used = count_trees_used + 1
        }

    }
    perc_ref = round(ref_count*100/tot_branches,2)
    depth =  round(mean(depths),2)
    return(list(ref_count,  perc_ref,count_trees_used,round(mean(relative_depths),2),depth))
}

method = "min" # because we want ALL to be above the threshold
me = compute_references(families, "MLsplits_MEbs", "", "", method, bs_threshold)
imd = compute_references(families, "MLsplits_IMDbs", "", "", method, bs_threshold)
ml = compute_references(families, "ML", "", "", method, bs_threshold)
imd_me = compute_references(families, "MLsplits_IMDbs", "MLsplits_MEbs", "", method, bs_threshold)
imd_ml = compute_references(families, "MLsplits_IMDbs", "ML", "", method, bs_threshold)
me_ml = compute_references(families, "MLsplits_MEbs", "ML", "", method, bs_threshold)
imd_me_ml = compute_references(families, "MLsplits_IMDbs", "MLsplits_MEbs", "ML", method, bs_threshold)


# S3
s3_colnames = c("Number ref. branches", "% branches", "Number of trees", "average relative depth", "average depth")
ref_table = rbind(imd,me, ml, imd_me, imd_ml, me_ml, imd_me_ml)
colnames(ref_table) = s3_colnames
rownames = gsub("_", "+", toupper(rownames(ref_table)))
#write table
ref_table = cbind(rownames, ref_table)
colnames(ref_table) = c("Ref branches", s3_colnames)
write.table(ref_table, paste(source_data, "../tables","S3.csv", sep = "/"), sep = ",", quote = F, row.names = F)





