setwd('~/coen_share/Hugo_analysis/parapatric_pools/popoolation/')

#################################
# Ploting Fst for a 50kb window #
#################################

#get data from saved tables
fst50 = read.table('./fst/magenta-yellow_ecori.window50.step10.fst')
rc = read.table('./allele_count/magenta-yellow_ecori.rc',na.strings="")

hyperplot = function(inputfst,inputrc,scaffold,region=FALSE,barcol=c("#FF8100","#0184A9"),linecol=c("black","#83D300","#D8006A")){
	
	meanfst = mean(inputfst[!is.na(inputfst[,6]),6])	#store the mean of whole-genome fst
	
	inputfst1 = inputfst[inputfst[,1]==scaffold,]	#from the inputfst select the ones refering to the scaffold of interest
	inputfst2 = inputfst1
	if(region){
		inputfst2 = inputfst1[inputfst1[,2]>region[1] & inputfst1[,2]<region[2],]	#from the inputfst select the ones refering to the region of interest
	}
	inputfst2 = droplevels(inputfst2)
	
	#Create the Fst plot
	par(fig=c(0,1,0.46,1))
	plot(inputfst2[,2]/1000,inputfst2[,6],xlab="",ylab="Fst",xaxt="n",las=1,type="l", ylim = c(0,1), main = scaffold, font.main="2",cex.main = 0.75)
	axis(3)
	abline(h=meanfst,lty=2,col="grey70")
	segments(min(inputfst2[,2])/1000,1,((max(inputfst2[,2])^2/1000)/max(inputfst[,2]))+(min(inputfst2[,2])/1000),1,lwd=5, col = "black",lend =3)
	abline(v=72000)
	#Prepare SNP count data
	rcscaff = inputrc[inputrc[,1]==scaffold,]	#select SNPs for selected scaffold
	rcscaff = droplevels(rcscaff)
	rcscaff = rcscaff[order(rcscaff[,2]),]	#make sure the table is ordered by position number
	#Shift things according to which is more frequent in either population
	rcscaff[rcscaff$ma_allele1 == rcscaff$ma_allele2 & rcscaff$ma_count1/rcscaff$depth1 < rcscaff$ma_count2/rcscaff$depth2,10] = (rcscaff$depth1 - rcscaff$ma_count1)[rcscaff$ma_allele1 == rcscaff$ma_allele2 & rcscaff$ma_count1/rcscaff$depth1 < rcscaff$ma_count2/rcscaff$depth2]
	rcscaff[rcscaff$ma_allele1 == rcscaff$ma_allele2 & rcscaff$ma_count2/rcscaff$depth2 < rcscaff$ma_count1/rcscaff$depth1,12] = (rcscaff$depth2 - rcscaff$ma_count2)[rcscaff$ma_allele1 == rcscaff$ma_allele2 & rcscaff$ma_count2/rcscaff$depth2 < rcscaff$ma_count1/rcscaff$depth1]


	#Barplots
	detach()
	attach(rcscaff)
	
	par(fig=c(0,1,0.37,0.70),new=T)
	plot(1:nrow(rcscaff),rep(0,nrow(rcscaff)),ylim=c(0,1),type="n",xlab="",ylab="",tck=0,xaxt="n", yaxt="n",frame.plot=FALSE,bty="n")
	

	par(xpd=NA)
	segments(1:nrow(rcscaff),0,(nrow(rcscaff)*rcscaff[,2])/max(inputfst2[,2]),1,col=linecol[1],lwd=1)
	segments(which(ma_count2/depth2>0 & ma_count2/depth2<1 & ma_count1/depth1<1),0,((nrow(rcscaff)*rcscaff[ma_count2/depth2>0 & ma_count2/depth2<1 & ma_count1/depth1<1,2]))/max(inputfst2[,2]),1,col=linecol[2],lwd=1)	#shared
	segments(which(ma_count2/depth2==1 & ma_count1/depth1==1),0,(nrow(rcscaff)*rcscaff[ma_count2/depth2==1 & ma_count1/depth1==1,2])/max(inputfst2[,2]),1,col=linecol[3],lwd=1)	#diagnostic
	
	par(fig=c(0,1,0.30,0.71),new=T)
	plot(1:nrow(rcscaff),rep(0,nrow(rcscaff)),ylim=c(0,1),type="n",xlab="",ylab="",tck=0,xaxt="n", yaxt="n",frame.plot=FALSE,bty="n")
	segments(1:nrow(rcscaff),0,1:nrow(rcscaff),0.5,col=linecol[1],lwd=2)
	segments(which(ma_count2/depth2>0 & ma_count2/depth2<1 & ma_count1/depth1<1),0,which(ma_count2/depth2>0 & ma_count2/depth2<1 & ma_count1/depth1<1),0.5,col=linecol[2],lwd=2)	#shared
	segments(which(ma_count2/depth2==1 & ma_count1/depth1==1),0,which(ma_count2/depth2==1 & ma_count1/depth1==1),0.5,col=linecol[3],lwd=2)	#diagnostic
	text((nrow(rcscaff)+8),0.1,as.character(length(which(ma_count2/depth2==1 & ma_count1/depth1==1))),cex=0.7,col=linecol[3])
	text((nrow(rcscaff)+8),0.25,as.character(length(which(ma_count2/depth2>0 & ma_count2/depth2<1 & ma_count1/depth1<1))),cex=0.7,col=linecol[2])
  text((nrow(rcscaff)+8),0.4,as.character(nrow(rcscaff)-(length(which(ma_count2/depth2>0 & ma_count2/depth2<1 & ma_count1/depth1<1))+length(which(ma_count2/depth2==1 & ma_count1/depth1==1)))),cex=0.7,col=linecol[1])


}

hyperplot(inputfst=fst50,inputrc=rc,scaffold="scaffold91")
hyperplot(inputfst=fst50,inputrc=rc,scaffold="scaffold91",region=c(0,10000))

plot(fst50[fst50$scaffold=="scaffold91",2],fst50[fst50$scaffold=="scaffold91",6],type="l")






























