# Load necessary libraries
library(stringr)

# Load args from command line
args = commandArgs(trailingOnly=TRUE)
# First argument is input gtf file
genepred_path = args[1]
# Second argument is list of gene names
sublist_path = args[2]
# Third argument is output file name
output_path = args[3]

# Read in input genepred file and list of gene names
orifile = read.table(genepred_path, header = FALSE, sep = '\t')
gene_sublist = read.table(sublist_path)

# Reformat names in list of gene names
gene_sublist = gene_sublist$V1
sublist_names = c()
for (i in 1:length(gene_sublist)) {
  sublist_names[i] = toString(gene_sublist[i])
}
sublist_names = sublist_names[order(sublist_names)]

# Reformat names in genepred file to match list of gene names exactly
all_names = c()
for (i in 1:length(orifile$V1)) {
  genepred_name = toString(orifile$V1[i])
  nametoadd = substr(genepred_name, 1, nchar(genepred_name)-2)
  all_names = append(all_names, nametoadd)
}

# Go through list of gene names to find overlap between list and input genepred file
overlap = data.frame()
for (i in 1:length(sublist_names)) {
  idx = match(sublist_names[i], all_names)
  if (is.na(idx) == FALSE) {
    overlap = rbind(overlap, orifile[idx,])
  }
}

# Write output as genepred file
write.table(overlap, file=output_path, col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
