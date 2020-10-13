# Load necessary libraries
library(stringr)

# Load args from command line
args = commandArgs(trailingOnly=TRUE)
# First argument is input gtf file
gtf_path = args[1]
# Second argument is list of gene names
sublist_path = args[2]
# Third argument is output file name
output_path = args[3]

# Read in input gtf and list of gene names
gtf = read.table(gtf_path, header = FALSE, sep = '\t')
gtf_with_quotes = read.table(gtf_path, header = FALSE, sep = '\t', quote="")
gene_sublist = read.table(sublist_path)

# Reformat names in list of gene names
gene_sublist = gene_sublist$V1
sublist_names = c()
for (i in 1:length(gene_sublist)) {
  sublist_names[i] = toString(gene_sublist[i])
}

# Reformat names in gtf to match list of gene names exactly
all_names = c()
for (i in 1:length(gtf$V9)) {
  geneinfo = str_split(toString(gtf$V9[i])," ", simplify=TRUE)
  nametoadd = geneinfo[1,2]
  all_names = append(all_names, substr(nametoadd, 1, (nchar(nametoadd)-3)))
}


# Go through list of gene names to find overlap between list and input gtf
overlap_gtf = data.frame()
for (i in 1:length(sublist_names)) {
  idx = match(sublist_names[i], all_names)
  if (is.na(idx) == FALSE) {
    overlap_gtf = rbind(overlap_gtf, gtf_with_quotes[idx,])
  }
}

# Write output as gtf file
write.table(overlap_gtf, file=output_path, col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
