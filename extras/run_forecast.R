if(!"RegTree"%in%installed.packages()){
  install.packages("~/Code/ac-portal/lib/RegTree", repos = NULL)
}

# Load needed libraries
library(RegTree)
library(dateutils)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

print(args)

df_path <- args[[1]]
lib_path <- args[[2]]
target <- args[[3]]
regress <- as.logical(args[[4]])

out <- 

write.csv(dtout, "/tmp/latest.csv", row.names = FALSE)










