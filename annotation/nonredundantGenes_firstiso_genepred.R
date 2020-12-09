# Load necessary libraries
library(stringr)

# Load args from command line
args = commandArgs(trailingOnly=TRUE)
# First argument is input gtf file
genepred_path = args[1]
# Second argument is output file name
output_name = args[2]

## Read in original genepred file
ori_genepred = read.table(genepred_path, header = FALSE, sep = '\t')
### Remove redundant rows of file by keeping the first row corresponding to each gene name
new_genepred = data.frame()
currname = 0
for (i in 1:length(ori_genepred$V1)) {
  nametoadd = toString(ori_genepred$V1[i])
  nametoadd = substr(nametoadd, 1, nchar(nametoadd)-2)
  if (nametoadd != currname) {
    new_genepred = rbind(new_genepred, ori_genepred[i,])
    currname = nametoadd
  }
}

## Write output as genepred file
write.table(new_genepred, file=output_name, col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

