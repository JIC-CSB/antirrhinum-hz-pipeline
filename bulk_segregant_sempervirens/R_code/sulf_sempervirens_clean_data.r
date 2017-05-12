###########################
# Bulk segregant analysis #
# Sempervirens SULF       #
###########################

#
# Fst
#
#Remove "="
fst50[,6] = as.numeric(gsub("[1,2,3,4]:[1,2,3,4]=","",fst50[,6]))


#give friendly collumn names to each of them
mycolnames = c("chr","pos","nSNP","cov","depth","mgt.ora")
colnames(fst50) = mycolnames

#remove windows with low coverage
fst50 = fst50[which(fst50$cov>0.01),]

#add linkage group and position
fst50 = merge(fst50, scaf[,c(1,2,4)], by.x="chr", by.y="scaf", all.x=T)
fst50$lg = as.character(fst50$lg)
fst50$lg[which(is.na(fst50$lg))] = "not_mapped"
fst50$lg = as.factor(fst50$lg)

#order table by linkage group
fst50 = fst50[order(fst50$lg, fst50$pos.cM, fst50$chr, fst50$pos),]

