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


get_splits_as_pp_binary <- function(splits, bs_threshold = 0){
    # Read the split file and return a binary vector of the splits that have a bootstrap value greater than the threshold
    # for each split we have a 1 if the split is greater than the threshold and 0 otherwise
    splits$V4 = as.numeric(splits$V4)
    binary_splits = as.numeric(splits$V4 >= bs_threshold)
    return(binary_splits)
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

combine_bs_3 <- function(splits, splits_2,splits_3, method){
    # order the splits by the split id
    splits = splits[order(splits$V5),]
    splits_2 = splits_2[order(splits_2$V5),]
    splits_3 = splits_3[order(splits_3$V5),]
    # check if we have the same splits, if not return an error
    if (length(setdiff(splits$V5, splits_2$V5)) > 0){
        stop("The splits are not the same")
    }
    if (method == "arithmetic_average"){
        splits$V4 = as.numeric(splits$V4)
        splits_2$V4 = as.numeric(splits_2$V4)
        splits_3$V4 = as.numeric(splits_3$V4)
        splits$V4 = (splits$V4 + splits_2$V4 + splits_3$V4)/3
    }
    else if (method == "geometric_average"){
        splits$V4 = as.numeric(splits$V4)
        splits_2$V4 = as.numeric(splits_2$V4)
        splits_3$V4 = as.numeric(splits_3$V4)
        splits$V4 = (splits$V4 * splits_2$V4 * splits_3$V4)^(1/3)
    }else if (method == "max"){
        splits$V4 = as.numeric(splits$V4)
        splits_2$V4 = as.numeric(splits_2$V4)
        splits_3$V4 = as.numeric(splits_3$V4)
        splits$V4 = pmax(splits$V4, splits_2$V4, splits_3$V4)
    }else if (method == "min"){
        splits$V4 = as.numeric(splits$V4)
        splits_2$V4 = as.numeric(splits_2$V4)
        splits_3$V4 = as.numeric(splits_3$V4)
        splits$V4 = pmin(splits$V4, splits_2$V4,  splits_3$V4)
    }
    return(splits)
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

get_splits_dir_names <- function(tree_type, ref , ncols = 200){

    ncols = as.character(ncols)
    if(ncols == "200"){
        ncols = ""
    }
    splits_dir = c()
    bs_trees = strsplit(ref,split = "+", fixed = TRUE)[[1]]

    for (bs_tree in bs_trees){
        if(bs_tree == tree_type && ncols == ""){
            splits_dir = c(splits_dir, tree_type)
        }else{
            splits_prefix = paste(tree_type, "splits", sep = "")
            splits_suffix = paste(bs_tree, "bs",ncols, sep = "")
            splits_dir = c(splits_dir, paste(splits_prefix, splits_suffix, sep = "_"))
        }
    }
    return(splits_dir)
}


get_proven_positives <- function(fam, tree_type, ref, bs_threshold){

    # extract split files 
    splits_dir = get_splits_dir_names(tree_type,ref,200)

    # combine bs 
    method = "min"
    if(length(splits_dir) == 2){
        combined_splits = combine_bs(parse_splits(fam, splits_dir[1], output_dir), parse_splits(fam, splits_dir[2], output_dir), method)
    }else if(length(splits_dir) == 3){
        combined_splits = combine_bs_3(parse_splits(fam, splits_dir[1], output_dir), parse_splits(fam, splits_dir[2], output_dir), parse_splits(fam, splits_dir[3], output_dir), method)
    }else{
        combined_splits = parse_splits(fam, splits_dir[1], output_dir)
    }

    # get the proven positives as a binary vector
    proven_positives = get_splits_as_pp_binary(combined_splits, bs_threshold)

    return(proven_positives)
}


get_splits_bs <- function(fam, tree_type, bs_type, ncols, method = ""){
    splits_dir = get_splits_dir_names(tree_type,bs_type,ncols)
    if(length(splits_dir) == 2){
        combined_splits = combine_bs(parse_splits(fam, splits_dir[1], output_dir), parse_splits(fam, splits_dir[2], output_dir), method)
    }else if(length(splits_dir) == 3){
        combined_splits = combine_bs_3(parse_splits(fam, splits_dir[1], output_dir), parse_splits(fam, splits_dir[2], output_dir), parse_splits(fam, splits_dir[3], output_dir), method)
    }else{
        combined_splits = parse_splits(fam, splits_dir[1], output_dir)
    }
    return(as.numeric(combined_splits$V4))
}


calculate_auc <- function(true_positives, predicted_positives) {
    # sort the predicted positives in descending order (and keep the true positives in the same order)
    n <- length(true_positives)
    data <- data.frame(score = predicted_positives, label = true_positives)
    data <- data[order(-data$score), ]
    best_mcc <- -1
    best_threshold <- 0

    # Initialize variables
    auc <- 0.0
    prev_fpr <- 0.0
    prev_tpr <- 0.0
    tp <- 0
    fp <- 0
    positives <- 0
    negatives <- 0

    # Count the number of positives and negatives

    for (i in 1:n) {
        if (data[i,]$label == 1) {
            positives <- positives + 1
        } else {
            negatives <- negatives + 1
        }
    }

    i <- 1
    while (i <= n) {
        j <- i
        tie_fp <- 0
        tie_tp <- 0
        
        # Handle ties by aggregating them
        while (j <= n && data[j,]$score == data[i,]$score) {
            if (data[j,]$label == 1) {
                tie_tp <- tie_tp + 1
            } else {
                tie_fp <- tie_fp + 1
            }
            j <- j + 1
        }
        
        tp <- tp + tie_tp
        fp <- fp + tie_fp
        
        tpr <- tp / positives
        fpr <- fp / negatives
        
        if (fp == 0 && tp > 0) {
            prev_tpr <- tpr
        } else if (i > 1 || (fp > 0 && tp > 0)) {
            auc <- auc + (fpr - prev_fpr) * (tpr + prev_tpr) / 2.0
            prev_fpr <- fpr
            prev_tpr <- tpr
        }
        
        i <- j  # Move to the next group of scores
    }
  
  return(list(auc = auc, best_threshold = best_threshold, best_mcc = best_mcc))
}


get_best_mcc <- function(predicted_values, proven_positives, thresholds){
    best_threshold = 0
    best_mcc = -1
    auc_v = 0 
    prev_tpr = 0
    prev_fpr = 0
    for( threshold in thresholds){
        predicted_labels <- ifelse(predicted_values > threshold, 1, 0)
        predicted_labels = factor(predicted_labels, levels = c(0,1))

        confusion_matrix = table(predicted_labels, proven_positives)
        tp = confusion_matrix[2,2]
        tn = confusion_matrix[1,1]
        fp = confusion_matrix[1,2]
        fn = confusion_matrix[2,1]
        tpr = tp / sum(proven_positives)
        fpr = fp / sum(proven_positives == 0)
        # if they are none, set them to 0
        if(is.nan(tpr)){
            tpr = 0
        }
        if(is.nan(fpr)){
            fpr = 0
        }

        mcc = (tp*tn - fp*fn) / sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
        if(is.nan(mcc)){
            mcc = -1
        }
        if(mcc > best_mcc){
            best_mcc = mcc
            best_threshold = threshold
            best_sensitivity = confusion_matrix[2,2] / sum(proven_positives)
            best_specificity = confusion_matrix[1,1] / sum(proven_positives == 0)
        }
    }
    
    return(list(best_mcc = best_mcc, best_threshold = best_threshold, best_sensitivity = best_sensitivity, best_specificity = best_specificity))
}

