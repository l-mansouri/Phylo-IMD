library(phangorn)

aligners=c('mTMalign')
trimming=c('untrimmed')

# ! PIPELINE OUTPUT FOLDER, CHANGE ACCORDINGLY
output_folder = '/home/luisasantus/Desktop/crg_cluster/NF_draft/'
# ! SOURCE DATA FOLDER, CHANGE ACCORDINGLY
source_data = '/home/luisasantus/Desktop/crg_cluster/projects/Phylo-IMD/analysis/source_data/'
tables = paste( source_data, "../tables/", sep = "/" )

df=read.table (paste(source_data,'/titration_every_5_reference_branches_ME+3d+ML.txt', sep = "/"))
colnames(df)=c('ncol','avg1d', 'sd1d', 'avg3d', 'sd3d', 'avgML', 'sdML', 'x')

thresholds = c(0.6, 0.75, 0.9)
colnames = c("avg1d", "avg3d", "avgML")
# for each threshold and each column name, get the ncol value and store in df with the corresponding threshold as column and column name as row
df_thresholds = data.frame(matrix(ncol = length(thresholds), nrow = length(colnames)))
colnames(df_thresholds) = thresholds
rownames(df_thresholds) = colnames
for (threshold in thresholds){
    for (colname in colnames){
        diff =abs(df[colname] - threshold)
        # remove everything negative
        # remove everything from diff that is above the threshold
        df_thresholds[colname,as.character(threshold)] = df[which.min(diff[[colname]]),"ncol"]
    }
}

# rename rows to Seq-ME, IMD-ME, Seq-ML
rownames(df_thresholds) = c("Seq-ME", "IMD-ME", "Seq-ML")
# rename columns to percentage
colnames(df_thresholds) = c("60%", "75%", "90%")
# reorder rows
df_thresholds = df_thresholds[c("Seq-ME", "Seq-ML","IMD-ME"),]
write.table(df_thresholds, file=paste(tables, paste("Table2.csv", sep=''), sep = "/"), append=F,quote=F, sep = ",", eol = "\n", na = "NA", dec = ".",row.names =F, col.names =T)

