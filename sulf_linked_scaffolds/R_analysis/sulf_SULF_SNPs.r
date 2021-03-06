##########################
# SNP frequency analysis #
##########################

rm(list=ls())

## Preparing the data ##

setwd('//group-data/shared/Research-Groups/Enrico-Coen/wgs/Project_GEL_EnricoCoen_JIC.EC.20120425.01/R_analysis')

# read in .sync file that was generated by the Popoolation2 script mpileup2sync
snpfreq = read.table('//group-data/shared/Research-Groups/Enrico-Coen/wgs/Project_GEL_EnricoCoen_JIC.EC.20120425.01/pileup/SULF_sulf.no_dupl.sync_rc')

# use the fasta index of the reference sequence to get scaffold lengths
scafflength = read.table('//group-data/shared/Research-Groups/Enrico-Coen/ref/snapdragon.scaf.version1.0.fai')

## Function to process popoolation sync file and calculate SNP frequencies ##
parsePoPoolationRC = function(input, depth) {
	
	# First keep only the data where there are only two alleles
	input = input[input[, 7] != "rc",]
	input = input[input[, 4] == 2,]
	
	
	maPop1 = as.numeric(unlist(strsplit(as.character(input[, 10]), "\\/")))	# major allele count in population 1
	maPop2 = as.numeric(unlist(strsplit(as.character(input[, 11]), "\\/")))	# major allele count in population 2
	maAlleles = as.character(unlist(strsplit(as.character(input[, 8]), "")))	# what the major alleles are in each population

	input[, 14] = maPop1[seq(1,nrow(input)*2, 2)]	# this is the major allele count in pop 1
	input[, 15] = maPop1[seq(2,nrow(input)*2, 2)]	# this is the depth in pop1
	input[, 16] = maPop2[seq(1,nrow(input)*2, 2)]	# this is the major allele count in pop 2
	input[, 17] = maPop2[seq(2,nrow(input)*2, 2)]	# this is the depth in pop2
	input[, 18] = maAlleles[seq(1,nrow(input)*2, 2)]	# this is the major allele state in pop 1
	input[, 19] = maAlleles[seq(2,nrow(input)*2, 2)]	# this is the major allele state in pop 2
	
	# Clean anything with depth < than indicated
	input = input[input[, 15] >= depth & input[, 17] >= depth,]
	
	# just output some of the collumns
	input = input[,c(1:5,8,9,18:19,14:17)]
	colnames(input) = c(
		"chr", "pos", "rc", "allele_count", "allele_states",
		"major_alleles(maa)", "minor_alleles(mia)",
		"ma_allele1", "ma_allele2",
		"ma_count1", "depth1", "ma_count2", "depth2"
	)

	return(input)
}

# Calculate SNP frequencies
rc = parsePoPoolationRC(snpfreq, 10)

# Make new collumns with frequency
rc["freq1"] = rc$ma_count1 / rc$depth1
rc["freq2"] = rc$ma_count2 / rc$depth2

## Function to calculate allele frequencies across sliding window ##
count_window = function(inputdata, length_scaffolds, windowsize, stepsize) {

	scaffolds = levels(inputdata[, 1]) #list of scaffolds

	scaffold = c() # name of the scaffold
	position = c() # position of the Fst point: middle of the window
	nrFixedSNP = c() # how many fixed snps in that window
	nrTotalSNP = c() # how many total snps in window

    # Loop through all scaffolds	
	for(i in 1:length(scaffolds)){
	
		scaftab = inputdata[inputdata[, 1] == scaffolds[i], ]	# table containing only values for the current scaffold
		scaflength = length_scaffolds[length_scaffolds == scaffolds[i], 2]	# stores the length of this scaffold
		cat(scaffolds[i], " (", i, "/", length(scaffolds), ") \n")	# prints some text to inform which scaffold is being processed

		# Make sure this window fits in the scaffold
        if (scaflength < windowsize) { 
		
			windowtab=scaftab 
			scaffold = c(scaffold, scaffolds[i])
			position = c(position, round(scaflength / 2))
			nrFixedSNP = c(nrFixedSNP, nrow(windowtab[which(windowtab$ma_allele2 != windowtab$rc & windowtab$freq2 == 1), ]))
			nrTotalSNP = c(nrTotalSNP, nrow(windowtab))
		
        # If the window doesn't fit truncate it at the end of the scaffold before making frequency calculation
        } else {

			for (startpos in seq(1, scaflength, stepsize)) {

				endpos = startpos+windowsize - 1
				
				if (round(mean(c(startpos, endpos))) < scaflength) { # the sliding window goes as long as its center is still in the scaffold.
				
					windowtab = scaftab[scaftab[, 2] >= startpos & scaftab[, 2] < endpos, ]
					scaffold = c(scaffold, scaffolds[i])
					position = c(position, round(mean(c(startpos, endpos))))
					nrFixedSNP = c(nrFixedSNP, nrow(windowtab[which(windowtab$ma_allele2!=windowtab$rc & windowtab$freq2==1),]))
					nrTotalSNP = c(nrTotalSNP, nrow(windowtab))
				}
		
			}
		
		}
	
	}
	
	output = data.frame(scaffold, position, nrFixedSNP,nrTotalSNP)
	colnames(output) = c("chr", "pos", "FixedSNP", "TotalSNP") 
	
	return(output)
}


## Applying different filters of SNP proportions ##

snpfilter = rc[which(rc$ma_allele2!=rc$rc), ]	# 209,256
snpfilter = snpfilter[which(snpfilter$freq2 == 1), ]	# 20,743
snpfilter = droplevels(snpfilter)

levels(rc$chr)	# ~5000 scaffolds

window1 = count_window(inputdata=rc, length_scaffolds=scafflength, windowsize=100000, stepsize=100000)
window2 = count_window(inputdata=rc, length_scaffolds=scafflength, windowsize=50000, stepsize=50000)
window3 = count_window(inputdata=rc, length_scaffolds=scafflength, windowsize=10000, stepsize=10000)

window1["ratio"] = window1$FixedSNP/window1$TotalSNP
window1[window1$TotalSNP>50 & window1$FixedSNP>0 & window1$ratio>0.8,]

scafflength[scafflength$V1=="scaffold2654",2]

#window1b = window1[which(window1$ratio!="NaN" & window1$TotalSNP>),
head(window1b[rev(order(window1b$ratio)), ], 40)

save.image(file='//group-data/shared/Research-Groups/Enrico-Coen/wgs/Project_GEL_EnricoCoen_JIC.EC.20120425.01/R_analysis/SNP_counts_Sulf.RData')

# check for allele bias towards reference (Stock 7)
rc_SULF_s7 = rc[rc$rc == rc$ma_allele1, ]
rc_sulf_s7 = rc[rc$rc == rc$ma_allele2, ]
rc_SULF_90 = rc[rc$freq1 > 0.88 & rc$freq1 < 0.92, ]
rc_sulf_90 = rc[rc$freq2 > 0.88 & rc$freq2 < 0.92, ]
rc_SULF_s7_90 = rc_SULF_s7[rc_SULF_s7$freq1 > 0.88 & rc_SULF_s7$freq1 < 0.92, ]
rc_sulf_s7_90 = rc_sulf_s7[rc_sulf_s7$freq2 > 0.88 &ck rc_sulf_s7$freq2 < 0.92, ]
slices <- c(nrow(rc_SULF_s7_90),  nrow(rc_SULF_90) -  nrow(rc_SULF_s7_90))
lbls <- c('stock 7 allele', 'other allele')
pie(slices, labels = lbls, main='SULF allele bias at 0.88 < freq < 0.92')
slices <- c(nrow(rc_sulf_s7_90),  nrow(rc_sulf_90) -  nrow(rc_sulf_s7_90))
pie(slices, labels = lbls, main='sulf allele bias at 0.88 < freq < 0.92')

window_SULF_s7 = count_window(inputdata=rc_SULF_s7, length_scaffolds=scafflength, windowsize=100000, stepsize=100000)
window_SULF_s7["ratio"] = window_SULF_s7$FixedSNP / window_SULF_s7$TotalSNP
window_sulf_s7 = count_window(inputdata=rc_sulf_s7, length_scaffolds=scafflength, windowsize=100000, stepsize=100000)

barplot(window1[ window1$TotalSNP > 100, ]$FixedSNP / window1[ window1$TotalSNP > 100, ]$TotalSNP)
par(new=T)
barplot(window_SULF_s7[ window_SULF_s7$TotalSNP > 100 & window_SULF_s7$FixedSNP > 0, ]$FixedSNP / window_SULF_s7[ window_SULF_s7$TotalSNP > 100 & window_SULF_s7$FixedSNP > 0, ]$TotalSNP, col = "royalblue")

rc_SULF_s7_intervals <- numeric(length = 50)
rc_sulf_s7_intervals <- numeric(length = 50)
intervals <- numeric(length = 50)
interval <- 0.02;

for (i in 0:49) {
	intervals[i] <- interval * i
	lower = i * interval
	upper = (i + 1) * interval
	rc_SULF_interval = rc[rc$freq1 > lower & rc$freq1 < upper, ]
	rc_sulf_interval = rc[rc$freq2 > lower & rc$freq2 < upper, ]
	rc_SULF_s7_interval <- rc_SULF_s7[rc_SULF_s7$freq1 > lower & rc_SULF_s7$freq1 < upper, ]
	rc_sulf_s7_interval <- rc_sulf_s7[rc_sulf_s7$freq2 > lower & rc_sulf_s7$freq2 < upper, ]
	rc_SULF_s7_intervals[i] = nrow(rc_SULF_s7_interval) / nrow(rc_SULF_interval) * 100 
	rc_sulf_s7_intervals[i] = nrow(rc_sulf_s7_interval) / nrow(rc_sulf_interval) * 100 
}

chart <- cbind(intervals, rc_SULF_s7_intervals)
barplot(chart[, 2], ylab="%age of stock7 as major allele in sulf pool", ylim=c(0,100), xlab="Allele Frequency", names.arg=chart[, 1], las=2)
