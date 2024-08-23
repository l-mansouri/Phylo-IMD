aligners=c('mTMalign')
trimming=c('untrimmed')
source_data = '/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data/'
plots = paste(source_data, '../plots/', sep = '')

fl=read.table(paste(source_data, 'list_of_families_with_all_rep_in_3d', sep = '/'))[,1]

folder_percid = '/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/results/mTMalign_percid/'

df = data.frame()
# iterate all files and extract line TOP 
for (fam in fl){
    print(fam)
    file = paste(folder_percid, fam, '.txt', sep = '')
    # read in txt file line by line 
    con <- file(file, "r")

    # Read the file line by line
    while(TRUE) {
    line <- readLines(con, n = 1, warn = FALSE)
    
    # Break the loop if we have reached the end of the file
    if (length(line) == 0) {
      break
    }

    # # Check if the line starts with "TOP"
    if (grepl("^TOP", line)) {
        # Split the line into columns based on whitespace
        cols <- strsplit(line, "\\s+")[[1]]

        # Extract the last three columns
        extracted_cols <- cols[(length(cols)-2):length(cols)]
        # add family name to the data frame
        extracted_cols = c(extracted_cols, fam)

        # Append the extracted columns to the data frame
        df <- rbind(df, as.data.frame(t(extracted_cols), stringsAsFactors = FALSE))
        
    }
    # add family name to the data frame
    
    }
    # close the file
    close(con)
}

# change column names
colnames(df) = c("seq1", "seq2", "percid", "family")
# make sure percid is numeric
df$percid = as.numeric(df$percid)
# save the data frame
write.table(df, paste(source_data, 'percids.txt', sep = ''), row.names = FALSE, quote = FALSE)
