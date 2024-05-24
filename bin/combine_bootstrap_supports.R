#!/usr/bin/env Rscript

library(phangorn)
library(argparser)

# Parse in arguments using flags using arg_parser
p <- arg_parser("This script computes the combined bootstrap support of two trees")
p <- add_argument(p, "--t1", help="The first tree file", type="character")
p <- add_argument(p, "--t2", help="The second tree file", type="character")
p <- add_argument(p, "--r1", help="The first replicates file", type="character")
p <- add_argument(p, "--r2", help="The second replicates file", type="character")
p <- add_argument(p, "--bs1", help="The first tree with bootstrap support", type="character")
p <- add_argument(p, "--bs2", help="The second tree with bootstrap support", type="character")
p <- add_argument(p, "--o", help="The combined tree with bootstrap support", type="character")
args <- parse_args(p)


# Read the trees
first_topology = read.tree(args$t1)
second_topology = read.tree(args$t2)
first_replicates = read.tree(args$r1)
second_replicates = read.tree(args$r2)


# Function to compute bootstrap support and add it to the tree file
compute_bs_support <- function(topology, replicates, output_file = ""){
  bs_support <- prop.clades(topology, replicates)
  bs_support[is.na(bs_support)]<-0
  topology$node.label = bs_support
  topology$node.label[1]=''
  if(output_file != ""){
    write.tree(topology, file = output_file)
  }
  return(topology)
}

# Compute first bootstrap support 
bs_first <- compute_bs_support(first_topology, first_replicates, args$bs1)
# Compute second bootstrap support 
bs_second <- compute_bs_support(second_topology, second_replicates, args$bs2)

#-------------------------------------------------------------------------------
#                   COMPUTE MULTISTRAP
#-------------------------------------------------------------------------------

second_replicates_on_first_topology <- compute_bs_support(first_topology, second_replicates)
external_support <- as.numeric(second_replicates_on_first_topology$node.label[-1])
first_support <- as.numeric(bs_first$node.label[-1])

# Compute the combination of the support values
supports <- cbind(first_support, external_support)
combined_supports <- c(c(""),  rowMeans(supports))

# Add the combined bootstrap support values to the first topology
first_topology$node.label <- combined_supports
write.tree(first_topology, file = args$o)




