###########################
# Bulk segregant analysis #
# Sempervirens SULF       #
###########################

#set working directory in windows
setwd('//group-data/shared/Research-Groups/Enrico-Coen/hybrid_zone/bulk_segregant/popoolation_outputs')

#set working directory in cluster
#setwd('/nbi/group-data/ifs/shared/Research-Groups/Enrico-Coen/hybrid_zone/bulk_segregant/popoolation_outputs/')

#source necessary scripts
source('//group-data/shared/Research-Groups/Enrico-Coen/sulf_mapping/bulk_segregant_sempervirens/R_code/sulf_sempervirens_load_data.r')
source('//group-data/shared/Research-Groups/Enrico-Coen/sulf_mapping/bulk_segregant_sempervirens/R_code/sulf_sempervirens_clean_data.r')


#
# Fst
#
#distribution
hist(fst50$mgt.ora, 
	breaks=200, main="", xlab="Magenta vs Orange", col="black")
text(0.3, 2000, paste("Median = ", median(fst50$mgt.ora), "\n99% = ", quantile(fst50$mgt.ora, 0.99)))


#
# Plotting Fst for whole-genome
#
#Change the yellow in the colour pallette to orange - or just use a palette of grey and black
palette(c("black", "grey","black", "grey","black", "grey","black", "grey","tan4"))

#make x coordinates
x = 1:nrow(fst50)

plot(x, fst50$mgt.ora, 
	pch=20, col=fst50$lg, xlab="", xaxt="n", ylab="Fst", main="Sempervirens SULF pools")
points(x[which(fst50$chr=="scaffold678-C8923637")], fst50$mgt.ora[which(fst50$chr=="scaffold678-C8923637")], pch=20, col="magenta")
points(x[which(fst50$chr=="scaffold117reverse")], fst50$mgt.ora[which(fst50$chr=="scaffold117reverse")], pch=20, col="magenta")
points(x[which(fst50$chr=="scaffold91")], fst50$mgt.ora[which(fst50$chr=="scaffold91")], pch=20, col="orange")
points(x[which(fst50$chr=="scaffold460")], fst50$mgt.ora[which(fst50$chr=="scaffold460")], pch=20, col="green")
points(x[which(fst50$chr=="scaffold316")], fst50$mgt.ora[which(fst50$chr=="scaffold316")], pch=20, col="purple")
points(x[which(fst50$chr=="scaffold455")], fst50$mgt.ora[which(fst50$chr=="scaffold455")], pch=20, col="pink")
legend("topleft", c("ROS scafs", "SULF scaf91", "scaf460", "monster scaf", "areusidin synthase", "unmapped scafs"), pch=20, col=c("magenta", "orange", "green", "purple","pink","tan4"))

points(x[which(fst50$chr=="scaffold776")], fst50$mgt.ora[which(fst50$chr=="scaffold776")], pch=20, col="red")


#
# Correlation with Fst from HZ
#
#Load Fst50 from hz
fst50hz = read.table('//group-data/shared/Research-Groups/Enrico-Coen/hybrid_zone/hybrid_zone_pools/popoolation_outputs/fst/yp4_yp2_yp1_mp2_mp4_mp11.stampy.markdup.fst50kb')
for(i in 6:20){
		fst50hz[,i] = as.numeric(gsub("[1,2,3,4,5,6,7,8]:[1,2,3,4,5,6,7,8]=","",fst50hz[,i]))
}
mycolnames = c("chr","pos","nSNP","cov","depth","yp4.yp2","yp4.yp1","yp4.mp2","yp4.mp4","yp4.mp11","yp2.yp1","yp2.mp2","yp2.mp4","yp2.mp11","yp1.mp2","yp1.mp4","yp1.mp11","mp2.mp4","mp2.mp11","mp4.mp11")
colnames(fst50hz) = mycolnames
#exclude windows with low coverage (<10%)
fst50hz[which(fst50hz$cov<0.1),6:20] = NA
#make an id
fst50hz["id"] = paste0(fst50hz$chr,"-",fst50hz$pos)
fst50["id"] = paste0(fst50$chr,"-",fst50$pos)

#merge Fst from bulk segregant data and hz
mergefst = merge(fst50hz[,c("chr", "pos", "yp4.mp11", "id")], fst50[,c("chr", "pos", "mgt.ora", "mgt.ora", "id")])

#plot correlation
plot(mergefst$yp4.mp11, mergefst$mgt.ora, xlab="HZ Fst outer pools", ylab="Sempervirens Fst SULF pools")
points(mergefst$yp4.mp11[which(mergefst$chr=="scaffold678-C8923637")], mergefst$mgt.ora[which(mergefst$chr=="scaffold678-C8923637")], pch=20, col="magenta")
points(mergefst$yp4.mp11[which(mergefst$chr=="scaffold117reverse")], mergefst$mgt.ora[which(mergefst$chr=="scaffold117reverse")], pch=20, col="magenta")
points(mergefst$yp4.mp11[which(mergefst$chr=="scaffold91")], mergefst$mgt.ora[which(mergefst$chr=="scaffold91")], pch=20, col="orange")
points(mergefst$yp4.mp11[which(mergefst$chr=="scaffold460")], mergefst$mgt.ora[which(mergefst$chr=="scaffold460")], pch=20, col="green")
points(mergefst$yp4.mp11[which(mergefst$chr=="scaffold316")], mergefst$mgt.ora[which(mergefst$chr=="scaffold316")], pch=20, col="purple")
points(mergefst$yp4.mp11[which(mergefst$chr=="scaffold455")], mergefst$mgt.ora[which(mergefst$chr=="scaffold455")], pch=20, col="pink")
legend("topright", c("ROS scafs", "SULF scaf", "scaf460", "monster scaf", "areusidin synthase"), pch=20, col=c("magenta", "orange", "green", "purple","pink"))

#plot scaffolds in case of interest
plot(fst50hz$pos[which(fst50hz$chr == "scaffold776")], fst50hz$yp4.mp11[which(fst50hz$chr == "scaffold776")], type="o", ylim=c(0,0.5))
abline(h=quantile(fst50hz$yp4.mp11, na.rm=T, c(0.9, 0.99)))


#
# Pull out top windows
#
top.windows = fst50[order(-fst50$mgt.ora),]	#order them from highest to lowest Fst
top.windows = top.windows[which(top.windows$mgt.ora > 0.2),]	#take a cutoff of Fst based visually on genome plot
top.windows = droplevels(top.windows)

#See how many scaffolds were picked
table(top.windows$chr)
dim(table(top.windows$chr))

#Pick those that have more than one window
table(top.windows$chr)[which(table(top.windows$chr) > 1)]

#make a table with: scaffold name, nr windows > 0.2, max Fst bulk, max Fst HZ, lg, cM

































