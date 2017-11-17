##Script to generate .tsv file of coverage per site split by read direction and proper paired status
## Edited 15/11/17

library(dplyr)

# Adjust as needed
setwd("~/Documents/work/jic/sulf/sulf_transposon_mapping/data")


#
# Functions ----
#
#Function to combine separate allele counts for variable sites into single integer

# Function to split values by "," and add them up
sumSplit <- function(x){
  # Split the values - this returns a list
  x <- strsplit(as.character(x), split = ",")
  
  # Loop the list and sum the values
  x <- lapply(x, function(i) sum(as.numeric(i)))
  
  # Unlist to return a vector
  unlist(x)
}


#
# Read mutant data ----
#
# Read data. Do not import strings as factors
paired_mut<-read.table("mutant.proper_paired.txt", sep="\t", header=F, stringsAsFactors = FALSE)
unpaired_mut<-read.table("mutant.mate_unmapped.txt", sep="\t", header=F, stringsAsFactors = FALSE)

# Assign column names
colnames(paired_mut)<-c("chr","pos","ref","alt","pDP","pADF","pADR")
colnames(unpaired_mut)<-c("chr","pos","ref","alt","uDP","uADF","uADR")

# Merge the two tables
reimport_mut <- merge(paired_mut, unpaired_mut, by=c("chr", "pos"), all = TRUE)
reimport_mut <- mutate_at(reimport_mut, vars(pDP:uADR), funs(sumSplit))

#convert na to zero coverage values
reimport_mut[is.na(reimport_mut)] <- 0

#calculate total depth per site
reimport_mut <- mutate(reimport_mut, total = pADF+pADR+uADF+uADR)

#write .tsv output
mock_tsv_mut <- select(chr, pos, total, pADF, pADR, uADF, uADR)
write.table(mock_tsv_mut, "mutant_counts.tsv", sep="\t", quote =F, col.names=F, row.names=F)


#
# Read WT data
#
# Read data. Do not import strings as factors
paired_wt <- read.table("WT.proper_paired.txt", sep="\t", header=F, stringsAsFactors = FALSE)
unpaired_wt <- read.table("WT.mate_unmapped.txt", sep="\t", header=F, stringsAsFactors = FALSE)

# Assign column names
colnames(paired_wt)<-c("chr","pos","ref","alt","pDP","pADF","pADR")
colnames(unpaired_wt)<-c("chr","pos","ref","alt","uDP","uADF","uADR")

# Merge the two tables
reimport_wt <- merge(paired_wt, unpaired_wt, by=c("chr", "pos"), all = TRUE)
reimport_wt <- mutate_at(reimport_wt, vars(pDP:uADR), funs(sumSplit))

#convert na to zero coverage values
reimport_wt[is.na(reimport_wt)] <- 0

#calculate total depth per site
reimport_wt <- mutate(reimport_wt, total = pADF+pADR+uADF+uADR)

#write .tsv output
mock_tsv_wt <- select(chr, pos, total, pADF, pADR, uADF, uADR)
write.table(mock_tsv_wt, "WT_counts.tsv", sep="\t", quote =F, col.names=F, row.names=F)



