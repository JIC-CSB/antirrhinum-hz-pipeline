###########################
# Bulk segregant analysis #
# Sempervirens SULF       #
###########################

#set working directory in windows
setwd('//group-data/shared/Research-Groups/Enrico-Coen/sulf_mapping/bulk_segregant_sempervirens/popoolation_outputs')

#set working directory in cluster
#setwd('/nbi/group-data/ifs/shared/Research-Groups/Enrico-Coen/hybrid_zone/bulk_segregant/popoolation_outputs/')

#
# FST
#
fst50 = read.table('./fst/SULF+_sulf.stampy_subs_0.05_ros.no_dupl.fst50kb')

#
# Scaffolds in linkage group
#
scaf = read.table('//group-data/shared/Research-Groups/Enrico-Coen/ref/Scaffold_In_Linkage_Group_JICv1.txt', header=T, sep="\t")

