
setwd('//group-data/shared/Research-Groups/Enrico-Coen/wgs/clean_data/map/')

#These tables are for reads with 0 < depth < 1000
sulf = read.table('ab-V660_R25.pair.v1.0.tsv2')
wt = read.table('ab-V661_R25.pair.v1.0.tsv2')

dim(wt)	#432,415,339 rows
dim(sulf)	#432,465,183 rows

#Change the collumn names to something more intuitive
colnames(wt) = c("chr","pos","depth","pairedF","pairedR","unpairF","unpairR")
colnames(sulf) = c("chr","pos","depth","pairedF","pairedR","unpairF","unpairR")

#Filter dataset
#Eliminate sites with missing data
wtf = wt[which(!is.na(wt$depth)),] #432,415,338
sulff = sulf[which(!is.na(sulf$depth)),] #432,465,183

hist(wtf$depth,breaks=seq(0,1000,1))
hist(sulff$depth,breaks=seq(0,1000,1))

#And eliminate sites with no unpaired reads
wtf = wtf[which(wtf$unpairF>0 | wtf$unpairR>0),]	#84,756,772
sulff = sulff[which(sulff$unpairF>0 | sulff$unpairR>0),]	#81,295,981

#And eliminate sites with no paired reads
wtf = wtf[which(wtf$pairedF>0 | wtf$pairedR>0),]	#76,972,460
sulff = sulff[which(sulff$pairedF>0 | sulff$pairedR>0),]	#73,810,980

#Filtering for putative insertions
wtf2 = wtf[which((wtf$unpairF>2 & wtf$pairedF<2 & wtf$unpairR<2 & wtf$pairedR>2) | (wtf$unpairF<2 & wtf$pairedF>2 & wtf$unpairR>2 & wtf$pairedR<2)),]	#8,146,507
sulff2 = sulff[which((sulff$unpairF>2 & sulff$pairedF<2 & sulff$unpairR<2 & sulff$pairedR>2) | (sulff$unpairF<2 & sulff$pairedF>2 & sulff$unpairR>2 & sulff$pairedR<2)),]	#7,891,060

#Create an identifier for each position
wtf2[,8] = paste(as.character(wtf2$chr),"-",as.character(wtf2$pos))
sulff2[,8] = paste(as.character(sulff2$chr),"-",as.character(sulff2$pos))

#Exclude sites common between wt and sulf
commonSites = match(sulff2$V8,wtf2$V8)

sulfSpecific = sulff2[which(is.na(commonSites)),] #1,881,284

#Exclude sites near start of scaffold
sulfSpecific = sulfSpecific[which(sulfSpecific$pos>=100),]	#1,834,231

#Do a classifier for which strand the missing pair is in
sulfSpecific[which(sulfSpecific$unpairF>1),9] = "F"
sulfSpecific[which(sulfSpecific$unpairR>1),9] = "R"

#Find inversion
myrows = 1:nrow(sulfSpecific)
sulfSpecific[which(sulfSpecific[myrows,]$V9=="F" & sulfSpecific[myrows+1,]$V9=="R"),10] = sulfSpecific[(which(sulfSpecific[myrows,]$V9=="F" & sulfSpecific[myrows+1,]$V9=="R"))+1,]$pos - sulfSpecific[which(sulfSpecific[myrows,]$V9=="F" & sulfSpecific[myrows+1,]$V9=="R"),]$pos

#How many inversions?
length(which(!is.na(sulfSpecific$V10)))	#17,634

#How many inversions across scaffolds?
length(which(sulfSpecific$V10<0))	#2,305

hist(sulfSpecific$V10,breaks=100000,xlim=c(0,5000),col="pink",main="Distance between F to R")

#Save tables
write.csv(sulfSpecific,'/net/nbi-cfs4.nbicluster/ifs/shared/Research-Groups/Enrico-Coen/sulf_mapping/sulf_transposon_mapping/data/sulfSpecific.csv', row.names = F, quote = F)

write.csv(sulff2,'/net/nbi-cfs4.nbicluster/ifs/shared/Research-Groups/Enrico-Coen/sulf_mapping/sulf_transposon_mapping/data/sulf_filtered_2.csv', row.names = F, quote = F)

write.csv(wtf2,'/net/nbi-cfs4.nbicluster/ifs/shared/Research-Groups/Enrico-Coen/sulf_mapping/sulf_transposon_mapping/data/wt_filtered_2.csv', row.names = F, quote = F)


##############################
# How long is each inversion?

#A function to calculate stuff for each inversion point
inversions = function(input=sulfSpecific){
	
	#Define the variables
	Iscaff = c()
	Fpos = c()
	Rpos = c()
	Flen = c()
	Rlen = c()
	
	#Define a list of scaffolds
	scaffolds = levels(input[,1]) 
	
	#Loop through each scaffold
	for(ii in 1:length(scaffolds)){
	
		myscaffold = input[which(input[,1]==scaffolds[ii]),]
	
		#Loop through table
		for(i in 1:nrow(myscaffold)){
			if(!is.na(myscaffold[i,10])){
				Iscaff = c(Iscaff,myscaffold[i,1])
				Fpos = c(Fpos,myscaffold[i,2])
				Rpos = c(Rpos,myscaffold[i+1,2])
				Flen = c(Flen,length(which(myscaffold[i-1000:i,9]=="F")))
				Rlen = c(Rlen,length(which(myscaffold[i:i+1000,9]=="R")))		
			}
		
		}
		
		#Print out a progress indicator
		cat("Scaffold",ii,"of",length(scaffolds),"\n")
	}
	
	output = data.frame(Iscaff,Ipos,Flen,Rlen)
	colnames(output) = c("chr","pos","Flen","Rlen")
	return(output)

}

myinversions = inversions()








