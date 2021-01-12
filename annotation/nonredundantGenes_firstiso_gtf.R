# Load necessary libraries
library(stringr)

# Load args from command line
args = commandArgs(trailingOnly=TRUE)
# First argument is input gtf file
gtf_path = args[1]
# Second argument is output file name
output_name = args[2]

## Read in original gtf file
gtf_with_quotes = read.table(gtf_path, header = FALSE, sep = '\t', quote="")
ori_gtf = read.table(gtf_path, header = FALSE, sep = '\t')
### Remove redundant rows of file by keeping the first row corresponding to each gene name
new_gtf = data.frame()
currname = 0
for (i in 1:length(ori_gtf$V1)) {
  geneinfo = str_split(toString(ori_gtf$V9[i])," ", simplify=TRUE)
  nametoadd = geneinfo[1,2]
  nametoadd = substr(nametoadd, 1, nchar(nametoadd)-3)
  if (nametoadd != currname) {
    new_gtf = rbind(new_gtf, gtf_with_quotes[i,])
    currname = nametoadd
  }
}

## Write output as genepred file
write.table(new_gtf, file=output_name, col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

