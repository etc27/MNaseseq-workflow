# Load necessary libraries
library(stringr)

# Load args from command line
args = commandArgs(trailingOnly=TRUE)
# First argument is input gtf file
gtf_path = args[1]
# Second argument is output file name
output_name = args[2]

## Read in gtf
gtf_with_quotes = read.table(gtf_path, header = FALSE, sep = '\t', quote="")
gtf = read.table(gtf_path, header = FALSE, sep = '\t')
### Remove redundant rows of gtf by keeping the row with each gene label ending in ".1"
# (e.g. keep only the row called AT1G01040.1 and remove all other rows with AT1G01040)
new_gtf = data.frame()
for (i in 1:length(gtf$V9)) {
  geneinfo = str_split(toString(gtf$V9[i])," ", simplify=TRUE)
  nametoadd = geneinfo[1,2]
  if (substr(nametoadd, nchar(nametoadd)-1, nchar(nametoadd)-1) == "1") {
    new_gtf = rbind(new_gtf, gtf_with_quotes[i,])
  }
}

## Write output as gtf file
write.table(new_gtf, file=output_name, col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

