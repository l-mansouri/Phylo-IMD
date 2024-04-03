# --------------------------------------------------
#  Compute the references 
# --------------------------------------------------
#!/usr/bin/env Rscript
library('phangorn')
library('geiger')


# FAMILIES 
families =read.table(paste("/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data",'list_of_families_for_titration', sep = "/"))[,1]

# BOOTSTRAP THRESHOLD
bs_threshold = 80 

# OUTPUT DIRECTORY OF SPLIT FILES 
output_dir = "/home/luisasantus/Desktop/crg_cluster/newphylo/NF_TMalign/mTMalign/titration_every_5_bootstrap_200_columns/split_files"

# internal namings - if using the pipeline this should not be touched as it is the default


# --------------------------------------------------
#  Functions
# --------------------------------------------------
parse_splits <- function(fam,splits_dir, output_dir){
    split_file_path <- list.files(paste(output_dir, splits_dir, sep = "/"), pattern = paste(fam,"_*", sep = ""), full.names = TRUE)
    split_file = read.csv(split_file_path, header = F, colClasses = "character")
    return(split_file)
}

get_splits <- function(splits, bs_threshold = 0){
    # Read the split file and return the splits that have a bootstrap value greater than the threshold
    number_splits = nrow(splits)
    splits$V4 = as.numeric(splits$V4)
    list_splits = as.data.frame(splits[splits$V4 >= bs_threshold,])
    if (nrow(list_splits) == 0){
        return(NULL)
    }
    else if (nrow(list_splits) == number_splits && bs_threshold != 0){
       return(NULL)
    }
    return(list_splits$V5)
}

get_depths <- function(splits, bs_threshold = 0){
    number_splits = nrow(splits)
    splits$V4 = as.numeric(splits$V4)
    splits$V3 = as.numeric(splits$V3)
    list_splits = as.data.frame(splits[splits$V4 >= bs_threshold,])
    if (nrow(list_splits) == 0){
        return(NULL)
    }
    else if (nrow(list_splits) == number_splits && bs_threshold != 0){
       return(NULL)
    }
    return(list_splits$V3)
}


combine_bs <- function(splits, splits_2, method){
    # order the splits by the split id
    splits = splits[order(splits$V5),]
    splits_2 = splits_2[order(splits_2$V5),]
    # check if we have the same splits, if not return an error
    if (length(setdiff(splits$V5, splits_2$V5)) > 0){
        stop("The splits are not the same")
    }
    if (method == "arithmetic_average"){
        splits$V4 = as.numeric(splits$V4)
        splits_2$V4 = as.numeric(splits_2$V4)
        splits$V4 = (splits$V4 + splits_2$V4)/2
    }
    else if (method == "geometric_average"){
        splits$V4 = as.numeric(splits$V4)
        splits_2$V4 = as.numeric(splits_2$V4)
        splits$V4 = sqrt(splits$V4 * splits_2$V4)
    }else if (method == "max"){
        splits$V4 = as.numeric(splits$V4)
        splits_2$V4 = as.numeric(splits_2$V4)
        splits$V4 = pmax(splits$V4, splits_2$V4)
    }else if (method == "min"){
        splits$V4 = as.numeric(splits$V4)
        splits_2$V4 = as.numeric(splits_2$V4)
        splits$V4 = pmin(splits$V4, splits_2$V4) 
    }
    return(splits)
}


bs_support_methods = c("ME", "IMD+ME", "IMD")

#--------------------------------------------------
# --------------------------------------------------
#                     Main
# --------------------------------------------------
# --------------------------------------------------

# --------------------------------------------------
#  SUPPLEMENTARY TABLE 3 
# --------------------------------------------------
method = "min"

get_ref_table_row <- function(families, dir, dir_2, dir_3, method, bs_threshold){
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
        imd_ref  = get_splits(splits_combined_bs, bs_threshold)

        #-------------------------
        # 2. Get the infos
        #-------------------------
        if(length(imd_ref) > 0){
            ref_count = ref_count + length(imd_ref)
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


me = get_ref_table_row(families, "MLsplits_MEbs", "", "", method, bs_threshold)
imd = get_ref_table_row(families, "MLsplits_IMDbs", "", "", method, bs_threshold)
ml = get_ref_table_row(families, "ML", "", "", method, bs_threshold)
imd_me = get_ref_table_row(families, "MLsplits_IMDbs", "MLsplits_MEbs", "", method, bs_threshold)
imd_ml = get_ref_table_row(families, "MLsplits_IMDbs", "ML", "", method, bs_threshold)
me_ml = get_ref_table_row(families, "MLsplits_MEbs", "ML", "", method, bs_threshold)
imd_me_ml = get_ref_table_row(families, "MLsplits_IMDbs", "MLsplits_MEbs", "ML", method, bs_threshold)

method

s3_colnames = c("Number ref. branches", "% branches", "Number of trees", "average relative depth", "average depth")
ref_table = rbind(imd,me, ml, imd_me, imd_ml, me_ml, imd_me_ml)
colnames(ref_table) = s3_colnames
ref_table




