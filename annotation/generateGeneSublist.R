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

# Reformat names in gene sublist
gene_sublist = gene_sublist$V1
sublist_names = c()
for (i in 1:length(gene_sublist)) {
  sublist_names[i] = toString(gene_sublist[i])
}

### Reformat names in full gene list to match subset gene names exactly
all_names = c()
for (i in 1:length(gtf$V9)) {
  geneinfo = str_split(toString(gtf$V9[i])," ", simplify=TRUE)
  nametoadd = geneinfo[1,2]
  all_names = append(all_names, substr(nametoadd, 1, (nchar(nametoadd)-3)))
}


### Go through gene sublist to find overlap
overlap_gtf = data.frame()
for (i in 1:length(sublist_names)) {
  idx = match(sublist_names[i], all_names)
  if (is.na(idx) == FALSE) {
    overlap_gtf = rbind(overlap_gtf, gtf_with_quotes[idx,])
  }
}

## Print gtf
write.table(overlap_gtf, file=output_path, col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
