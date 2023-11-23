setwd("C:/Users/aki/Desktop/P.borealis/P.borealis/")
require(ggplot2)
require(vcfR)
require(cowplot)
require(ggplot2)
require(ggrepel)
require(seqinr)

#require(BiocManager)
#BiocManager::install("qvalue")
#require(devtools)
#devtools::install_github("whitlock/OutFLANK")
require(OutFLANK)

#BiocManager::install("SNPRelate")
require(SNPRelate)

#remotes::install_github("privefl/bigsnpr")
require(bigsnpr)

#install.packages("robust")
require(robust)

#install.packages("bigstatsr")
require(bigstatsr)

###############################################################
# Pre-filtration (population and SNP level) data visualization 
###############################################################
# x<-read.vcfR("P.borealis_stacks.vcf")
# 
# queryMETA(x)
# queryMETA(x, element = 'DP')
# 
# ad <- extract.gt(x, element = "AD", as.numeric=TRUE)
# 
# dp <- extract.gt(x, element = "DP", as.numeric=TRUE)
# dp2<-data.frame(colnames(dp),colMeans(dp,na.rm = T))
# colnames(dp2)<-c("Sample","MeanDP")
# dp2$Max <- apply(dp, 2, function(x) max(x, na.rm = TRUE))
# dp2$Min <- apply(dp, 2, function(x) min(x, na.rm = TRUE))
# dp2$Median <- apply(dp, 2, function(x) median(x, na.rm = TRUE))
# 
# ggplot(dp2)+
#   geom_col(aes(x=Sample,y=MeanDP))+
#   geom_point(aes(x=Sample,y=Median))+
#   #geom_pointrange(aes(x=Sample,y=MeanDP, ymin=Min,ymax=Max))+
#   ylab("DP")+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6.5))
# 
# dp3<-data.frame(row.names(dp),rowMeans(dp,na.rm = T))
# colnames(dp3)<-c("SNP","MeanDP")
# dp3$MedianDP <- apply(dp, 1, function(x) median(x, na.rm = TRUE))
# dp3$Max <- apply(dp, 1, function(x) max(x, na.rm = TRUE))
# dp3$Min <- apply(dp, 1, function(x) min(x, na.rm = TRUE))
# 
# ggplot(dp3)+
#   geom_col(aes(x=SNP,y=MeanDP))+
#   geom_point(aes(x=SNP,y=MedianDP))+
#   #geom_pointrange(aes(x=Sample,y=MeanDP, ymin=Min,ymax=Max))+
#   ylab("DP")+
#   ylim(0,50)+
#   theme_classic()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# mean((dp3$MeanDP))
# sd((dp3$MeanDP))
# 
# ggplot(dp3[dp3$MeanDP>14 & dp3$MeanDP<22,])+
#   geom_col(aes(x=SNP,y=MeanDP))+
#   geom_point(aes(x=SNP,y=MedianDP))+
#   #geom_pointrange(aes(x=Sample,y=MeanDP, ymin=Min,ymax=Max))+
#   ylab("DP")+
#   ylim(0,50)+
#   theme_classic()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# hist((dp3[dp3$MeanDP>14 & dp3$MeanDP<22,]$MeanDP))
# 
# GQ <- extract.gt(x, element = "GQ", as.numeric=TRUE)
# 
# 
# indmiss<-read.table("out.imiss",sep="\t",header = T)
# indmiss$INDV
# 
# ind<-c("SU18_10", "SU18_11", "SU18_13", "SU18_17", "SU18_18", "SU18_19", "SU18_20", "SU18_21", "SU18_22", "SU18_23", "SU18_24", "SU18_3", "SU18_6", "SU18_7", "SU18_9", "AR18_26", "AR18_29", "AR18_30", "AR18_31", "AR18_32", "AR18_33", "AR18_34", "AR18_36", "AR18_39", "AR18_41", "AR18_42", "AR18_43", "AR18_44", "AR18_45", "AR18_47", "AR18_48", "SI18_49", "SI18_51", "SI18_52", "SI18_53", "SI18_56", "SI18_57", "SI18_58", "SI18_59", "SI18_60", "SI18_66", "SI18_67", "SI18_68", "SI18_70", "SI18_71", "SI18_72", "SD18_74", "SD18_75", "SD18_76", "SD18_78", "SD18_79", "SD18_80", "SD18_81", "SD18_84", "SD18_86", "SD18_87", "SD18_88", "SD18_89", "SD18_90", "SD18_91", "SD18_92", "SD18_94", "YS21_1", "YS21_10", "YS21_11", "YS21_12", "YS21_2", "YS21_3", "YS21_4", "YS21_5", "YS21_6", "YS21_7", "YS21_8", "YS21_9", "IS21_66", "IS21_67", "IS21_68", "IS21_69", "IS21_70", "IS21_71", "IS21_72", "IS21_73", "IS21_74", "IS21_75", "IS21_76", "UH21_141", "UH21_142", "UH21_143", "UH21_144", "UH21_145", "UH21_146", "UH21_147", "UH21_148", "UH21_149", "UH21_150", "UH21_151")
# samp<-c("10_SU18","11_SU18","13_SU18","17_SU18","18_SU18","19_SU18","20_SU18","21_SU18","22_SU18","23_SU18","24_SU18","3_SU18","6_SU18","7_SU18","9_SU18","26_AR18","29_AR18","30_AR18","31_AR18","32_AR18","33_AR18","34_AR18","36_AR18","39_AR18","41_AR18","42_AR18","43_AR18","44_AR18","45_AR18","47_AR18","48_AR18","49_SI18","51_SI18","52_SI18","53_SI18","56_SI18","57_SI18","58_SI18","59_SI18","60_SI18","66_SI18","67_SI18","68_SI18","70_SI18","71_SI18","72_SI18","74_SD18","75_SD18","76_SD18","78_SD18","79_SD18","80_SD18","81_SD18","84_SD18","86_SD18","87_SD18","88_SD18","89_SD18","90_SD18","91_SD18","92_SD18","94_SD18","1_YS21","10_YS21","11_YS21","12_YS21","2_YS21","3_YS21","4_YS21","5_YS21","6_YS21","7_YS21","8_YS21","9_YS21","66_IS21","67_IS21","68_IS21","69_IS21","70_IS21","71_IS21","72_IS21","73_IS21","74_IS21","75_IS21","76_IS21","141_UH21","142_UH21","143_UH21","144_UH21","145_UH21","146_UH21","147_UH21","148_UH21","149_UH21","150_UH21","151_UH21")
# pop<-c("SU18", "SU18", "SU18", "SU18", "SU18", "SU18", "SU18", "SU18", "SU18", "SU18", "SU18", "SU18", "SU18", "SU18", "SU18", "AR18", "AR18", "AR18", "AR18", "AR18", "AR18", "AR18", "AR18", "AR18", "AR18", "AR18", "AR18", "AR18", "AR18", "AR18", "AR18", "SI18", "SI18", "SI18", "SI18", "SI18", "SI18", "SI18", "SI18", "SI18", "SI18", "SI18", "SI18", "SI18", "SI18", "SI18", "SD18", "SD18", "SD18", "SD18", "SD18", "SD18", "SD18", "SD18", "SD18", "SD18", "SD18", "SD18", "SD18", "SD18", "SD18", "SD18", "YS21", "YS21", "YS21", "YS21", "YS21", "YS21", "YS21", "YS21", "YS21", "YS21", "YS21", "YS21", "IS21", "IS21", "IS21", "IS21", "IS21", "IS21", "IS21", "IS21", "IS21", "IS21", "IS21", "UH21", "UH21", "UH21", "UH21", "UH21", "UH21", "UH21", "UH21", "UH21", "UH21", "UH21")
# year<-c(rep("2018",62),rep("2021",34))
# indmiss$Sample<-ind
# indmiss$Sample2<-samp
# indmiss$Pop<-pop
# indmiss$Year<-year
# indmiss<-indmiss[order(indmiss$F_MISS,decreasing = T),]
# indmiss$Sample2<-factor(indmiss$Sample2,levels = indmiss$Sample2)
# ggplot(data=indmiss)+
#   geom_col(aes(x = Sample2, y = (F_MISS),fill=Year))+
#   ylab("Fraction of missing SNPs")+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = -45, vjust = 0.1,hjust = 0.1))
# 
# aggregate(F_MISS~Pop,indmiss,FUN=mean)

#New filtered file
x<-read.vcfR("../P.borealis_stacks.g5mac3dp10ind50g95maf05mdp20p1.recode.vcf")
queryMETA(x)
dp <- extract.gt(x, element = "DP", as.numeric=TRUE)
GQ <- extract.gt(x, element = "GQ", as.numeric=TRUE)

summary(dp)
hist(dp)
min(apply(dp, 2, function(x) min(x, na.rm = TRUE)))
max(apply(dp, 2, function(x) max(x, na.rm = TRUE)))

summary(GQ)
hist(GQ)
min(apply(GQ, 2, function(x) min(x, na.rm = TRUE)))
max(apply(GQ, 2, function(x) max(x, na.rm = TRUE)))


geno <- extract.gt(x) # Character matrix Containing the genotypes
position <- getPOS(x) # Positions in bp
chromosome <- getCHROM(x) # Chromosome information
pos_loc<-paste(chromosome,position,sep="_")

colnames(geno)
table(as.vector(geno[,1:9])) #SU18
table(as.vector(geno[,10:19])) #AR18
table(as.vector(geno[,20:29])) #SI18
table(as.vector(geno[,30:36])) #SD18
table(as.vector(geno[,37:47])) #YS21
table(as.vector(geno[,48:57])) #IS21
table(as.vector(geno[,58:68])) #UH21

pop2<-c(rep("SU18",9),rep("AR18",10),rep("SI18",10),rep("SD18",7),rep("YS21",11),rep("IS21",10),rep("UH21",11))
samp2<-c("11_SU18","18_SU18","19_SU18","20_SU18","21_SU18","22_SU18","24_SU18","3_SU18","9_SU18","29_AR18","30_AR18","31_AR18","36_AR18","39_AR18","42_AR18","43_AR18","44_AR18","45_AR18","48_AR18","49_SI18","51_SI18","52_SI18","57_SI18","59_SI18","66_SI18","67_SI18","68_SI18","70_SI18","72_SI18","74_SD18","75_SD18","80_SD18","81_SD18","88_SD18","89_SD18","92_SD18","1_YS21","10_YS21","11_YS21","12_YS21","2_YS21","3_YS21","4_YS21","5_YS21","6_YS21","8_YS21","9_YS21","67_IS21","68_IS21","69_IS21","70_IS21","71_IS21","72_IS21","73_IS21","74_IS21","75_IS21","76_IS21","141_UH21","142_UH21","143_UH21","144_UH21","145_UH21","146_UH21","147_UH21","148_UH21","149_UH21","150_UH21","151_UH21")

#prepare sample for OutFLANK analysis
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno),dimnames = list(pos_loc,colnames(geno)))

G[geno %in% c("0/0")] <- 0
G[geno %in% c("0/1", "1/0")] <- 1
G[geno %in% c("1/1")] <- 2

sum(is.na(geno))/(nrow(geno)*ncol(geno))

G[is.na(G)]<-9

#separate out the temporal samples
geno_18<-geno[,c(1:36)]
colnames(geno_18)

geno_21<-geno[,c(37:68)]
colnames(geno_21)

sum(is.na(geno_18))/(nrow(geno_18)*ncol(geno_18))

sum(is.na(geno_21))/(nrow(geno_21)*ncol(geno_21))

nll_18<-NULL
for(i in 1:nrow(geno_18)){
  if(length(unique(geno_18[i,]))==1){
    nll_18<-c(nll_18,i)
  }
}

nll_21<-NULL
for(i in 1:nrow(geno_21)){
  if(length(unique(geno_21[i,]))==1){
    nll_21<-c(nll_21,i)
  }
}

G18 <- matrix(NA, nrow = nrow(geno_18), ncol = ncol(geno_18),dimnames = list(pos_loc,colnames(geno_18)) )

G18[geno_18 %in% c("0/0")] <- 0
G18[geno_18 %in% c("0/1", "1/0")] <- 1
G18[geno_18 %in% c("1/1")] <- 2

G18[is.na(G18)]<-9

G21 <- matrix(NA, nrow = nrow(geno_21), ncol = ncol(geno_21),dimnames = list(pos_loc,colnames(geno_21)) )

G21[geno_21 %in% c("0/0")] <- 0
G21[geno_21 %in% c("0/1", "1/0")] <- 1
G21[geno_21 %in% c("1/1")] <- 2

G21[is.na(G21)]<-9

SNPmat18<-(t(G18))

SNPmat21<-(t(G21))


#modified function MakeDiploidFSTMat
MakeDiploidFSTMat<-function(SNPmat,locusNames,popNames){
  locusname <- unlist(locusNames)
  popname <- unlist(popNames)
  snplevs <- levels(as.factor(unlist(SNPmat)))
  if(any(!(snplevs%in%c(0,1,2,9)))==TRUE) {
    print("Error: Your snp matrix has a character other than 0,1,2 or 9")
    break
  }
  if (dim(SNPmat)[1] != length(popname)) {
    print("Error: your population names do not match your SNP matrix")
    break
  }
  if (dim(SNPmat)[2] != length(locusname)) {
    print("Error:  your locus names do not match your SNP matrix")
    break
  }
  writeLines("Calculating FSTs, may take a few minutes...")
  nloci <- length(locusname)
  FSTmat <- matrix(NA, nrow = nloci, ncol = 8)
  for (i in 1:nloci) {
    FSTmat[i, ] = unlist(getFSTs_diploids(popname, SNPmat[,i]))
    if (i%%10000 == 0) {
      print(paste(i, "done of", nloci))
    }
  }
  outTemp = as.data.frame(FSTmat)
  outTemp = cbind(locusname, outTemp)
  colnames(outTemp) = c("LocusName", "He", "FST", "T1", "T2", 
                        "FSTNoCorr", "T1NoCorr", "T2NoCorr", "meanAlleleFreq")
  return(outTemp)
}

getFSTs_diploids = function(popNameList, SNPDataColumn){  
  #eliminating the missing data for this locus
  popnames=unlist(as.character(popNameList))
  popNameTemp=popnames[which(SNPDataColumn!=9)]
  snpDataTemp=SNPDataColumn[SNPDataColumn!=9]
  
  HetCounts <- tapply(snpDataTemp, list(popNameTemp,snpDataTemp), length)
  HetCounts[is.na(HetCounts)] = 0
  
  #Case: all individuals are genetically identical at this locus
  if(dim(HetCounts)[2]==1){
    return (list(He=NA,FST=NA, T1=NA, T2=NA,FSTNoCorr=NA, T1NoCorr=NA, T2NoCorr=NA,meanAlleleFreq = NA))
  }
  
  if(dim(HetCounts)[2]==2){
    if(paste(colnames(HetCounts),collapse="")=="01"){HetCounts=cbind(HetCounts,"2"=0)}
    if(paste(colnames(HetCounts),collapse="")=="12"){HetCounts=cbind("0"=0,HetCounts)} 
    if(paste(colnames(HetCounts),collapse="")=="02"){HetCounts=cbind(HetCounts[,1],"1"=0, HetCounts[,2])}
  }
  
  out = WC_FST_Diploids_2Alleles(HetCounts)	
  return(out)
}

pop_18 <- pop2[c(1:36)]
pop_21 <- pop2[c(37:68)]

my_fst <- MakeDiploidFSTMat(t(G), locusNames = pos_loc, popNames = pop2)
my_fst_18 <- MakeDiploidFSTMat(t(G18), locusNames = pos_loc, popNames = pop_18)
my_fst_21 <- MakeDiploidFSTMat(t(G21), locusNames = pos_loc, popNames = pop_21)

plot(my_fst$He, my_fst$FST, xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))
plot(my_fst_18$He, my_fst_18$FST, xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))
plot(my_fst_21$He, my_fst_21$FST, xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))

# If a locus has a much lower sample size compared to the rest, it could have a broader error distribution (and therefore incorrectly inferred as an outlier).
plot(my_fst$FST, my_fst$FSTNoCorr)
abline(a=0,b=1,col="red")

plot(my_fst_18$FST, my_fst_18$FSTNoCorr)
abline(a=0,b=1,col="red")

plot(my_fst_21$FST, my_fst_21$FSTNoCorr)
abline(a=0,b=1,col="red")

out_trim <- OutFLANK(my_fst, NumberOfSamples=length(unique(pop2)), qthreshold = 0.001, Hmin = 0.1)
out_trim_18 <- OutFLANK(my_fst_18, NumberOfSamples=length(unique(pop_18)), qthreshold = 0.001, Hmin = 0.1)
out_trim_21 <- OutFLANK(my_fst_21, NumberOfSamples=length(unique(pop_21)), qthreshold = 0.001, Hmin = 0.1)

str(out_trim)
head(out_trim$results)
summary(out_trim$results$OutlierFlag)
summary(out_trim$results$pvaluesRightTail)

OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
## Zoom in on right tail
OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

hist(out_trim$results$pvaluesRightTail)

P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)

my_out <- which(P1$OutlierFlag==TRUE)
length(my_out)
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1),xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))
points(P1$He[my_out], P1$FST[my_out], col="blue",pch=20)

plot(as.factor(P1[!is.na(P1$He),]$LocusName[P1[!is.na(P1$He),]$He>0.1]), P1[!is.na(P1$He),]$FST[P1[!is.na(P1$He),]$He>0.1],
     xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
points(as.factor(P1[!is.na(P1$He),]$LocusName[P1[!is.na(P1$He),]$He>0.1])[my_out], P1[!is.na(P1$He),]$FST[P1[!is.na(P1$He),]$He>0.1][my_out], col="magenta", pch=20)  

OutLoc<-P1[P1$OutlierFlag==TRUE,]
OutLoc<-OutLoc[!is.na(OutLoc$OutlierFlag),]

str(out_trim_18)
head(out_trim_18$results)
summary(out_trim_18$results$OutlierFlag)
summary(out_trim_18$results$pvaluesRightTail)

OutFLANKResultsPlotter(out_trim_18, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
## Zoom in on right tail
OutFLANKResultsPlotter(out_trim_18 , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

hist(out_trim_18$results$pvaluesRightTail)

P1_18 <- pOutlierFinderChiSqNoCorr(my_fst_18, Fstbar = out_trim_18$FSTNoCorrbar, dfInferred = out_trim_18$dfInferred, qthreshold = 0.05, Hmin=0.1)

my_out_18 <- which(P1_18$OutlierFlag==TRUE)
length(my_out_18)
plot(P1_18$He, P1_18$FST, pch=19, col=rgb(0,0,0,0.1),xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))
points(P1_18$He[my_out_18], P1_18$FST[my_out_18], col="blue",pch=20)

plot(as.factor(P1_18[!is.na(P1_18$He),]$LocusName[P1_18[!is.na(P1_18$He),]$He>0.1]), P1_18[!is.na(P1_18$He),]$FST[P1_18[!is.na(P1_18$He),]$He>0.1],
     xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
points(as.factor(P1_18[!is.na(P1_18$He),]$LocusName[P1_18[!is.na(P1_18$He),]$He>0.1])[my_out_18], P1_18[!is.na(P1_18$He),]$FST[P1_18[!is.na(P1_18$He),]$He>0.1][my_out_18], col="magenta", pch=20)  

OutLoc18<-P1_18[P1_18$OutlierFlag==TRUE,]
OutLoc18<-OutLoc18[!is.na(OutLoc18$OutlierFlag),]


str(out_trim_21)
head(out_trim_21$results)
summary(out_trim_21$results$OutlierFlag)
summary(out_trim_21$results$pvaluesRightTail)

OutFLANKResultsPlotter(out_trim_21, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
## Zoom in on right tail
OutFLANKResultsPlotter(out_trim_21 , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

hist(out_trim_21$results$pvaluesRightTail)

P1_21 <- pOutlierFinderChiSqNoCorr(my_fst_21, Fstbar = out_trim_21$FSTNoCorrbar, dfInferred = out_trim_21$dfInferred, qthreshold = 0.05, Hmin=0.1)

my_out_21 <- P1_21$OutlierFlag==TRUE
sum(my_out_21,na.rm = T)
plot(P1_21$He, P1_21$FST, pch=19, col=rgb(0,0,0,0.1),xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))
points(P1_21$He[my_out_21], P1_21$FST[my_out_21], col="blue")

plot(as.factor(P1_21[!is.na(P1_21$He),]$LocusName[P1_21[!is.na(P1_21$He),]$He>0.1]), P1_21[!is.na(P1_21$He),]$FST[P1_21[!is.na(P1_21$He),]$He>0.1],
     xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
points(as.factor(P1_21[!is.na(P1_21$He),]$LocusName[P1_21[!is.na(P1_21$He),]$He>0.1])[my_out_21], P1_21[!is.na(P1_21$He),]$FST[P1_21[!is.na(P1_21$He),]$He>0.1][my_out_21], col="magenta", pch=20)  

OutLoc21<-P1_21[P1_21$OutlierFlag==TRUE,]
OutLoc21<-OutLoc21[!is.na(OutLoc21$OutlierFlag),]


OutLoc
OutLoc18
OutLoc21
#PCA
#vcf.fn<-"P.borealis_stacks.g5mac3dp10ind50g95maf05mdp20p1.recode.vcf"
#snpgdsVCF2GDS(vcf.fn, "P.borealis.gds",method = "copy.num.of.ref")

Pbor<-snpgdsOpen("../P.borealis.gds")
#snpgdsClose(Slam)

pcaC <- snpgdsPCA(Pbor,autosome.only=F)

pc.percent <- pcaC$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pcaC$sample.id,
                  EV1 = pcaC$eigenvect[,1],
                  EV2 = pcaC$eigenvect[,2],
                  EV3 = pcaC$eigenvect[,3],
                  EV4 = pcaC$eigenvect[,4],
                  #EV5 = pcaC$eigenvect[,5],
                  #EV6 = pcaC$eigenvect[,6],
                  stringsAsFactors = FALSE)
head(tab)

Per_exp<-head(round(pc.percent, 2))

year<-c(rep("2018",36),rep("2021",32))
tab$pop<-pop2
tab$year<-year
tab$ind<-samp2 #samp2 object from line 138
cl<-c("blue","orange","cyan", "lightgreen", "darkgreen","red", "black")
#cl<-c("blue","cyan", "orange", "darkgreen","red")
#okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E15759")
shp<-c(0,1,2,3,4,5,6)
#shp<-c(0,1,2,3,4)
ggplot(tab, aes(EV1,EV2,color=pop,shape=pop)) +
  xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) + 
  ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + 
  stat_ellipse(level=0.75,size=2)+
  geom_point(size=5) +
  #  scale_shape_manual(name="Pop", labels=unique(pop)[order(unique(pop))], values=shp)+
  #  scale_color_manual(name="Pop", labels=unique(pop)[order(unique(pop))], values=cl)+
  scale_shape_manual(name="Pop", labels=unique(pop2)[order(unique(pop2))], values=shp)+
  scale_color_manual(name="Pop", labels=unique(pop2)[order(unique(pop2))], values=cl)+
  labs(color="") + 
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

ggplot(tab, aes(EV1,EV2,color=year,shape=year)) +
  xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) + 
  ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + 
  geom_point(size=3,stroke=1.2) + 
  stat_ellipse(level=0.75,size=1)+
  scale_color_manual(name="Year", labels=unique(year)[order(unique(year))],values = c("red","cyan"))+
  scale_shape_manual(name="Year", labels=unique(year)[order(unique(year))], values= c(0,6))+
  labs(color="") + 
  theme_classic()

#smartPCA
spca_tab<-read.table("P.borealis_stacks.g5mac3dp10ind50g95maf05mdp20p1.recode.evec")
spca_Per_exp<-read.table("P.borealis_stacks.g5mac3dp10ind50g95maf05mdp20p1.recode.eval")
spca_tab$V7<-pop2
colnames(spca_tab)<-c("sample.id","EV1","EV2","EV3", "EV4","EV5","Pop")
#spca_tab$order<- c(rep(1,25),rep(2,12),rep(3,5),rep(4,11),rep(9,3),rep(5,12),rep(6,8),rep(7,2),rep(8,13),10)
#spca_tab<-spca_tab[order(spca_tab$order),]

barplot(spca_Per_exp$V1)

ggplot(spca_tab,aes(EV1,EV2,color=Pop,shape=Pop)) +
  xlab(paste("PC1 (",round(spca_Per_exp[1,],2),"%)",sep="")) + 
  ylab(paste("PC2 (",round(spca_Per_exp[2,],2),"%)",sep="")) + 
  stat_ellipse(level=0.75,size=2)+
  geom_point(size=5) +
  #geom_point(size=3, subset = .(label == 'point'))+
  scale_shape_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=shp)+
  scale_color_manual(name="Pop", labels=unique(spca_tab$Pop)[order(unique(spca_tab$Pop))], values=cl)+
  labs(color="") + 
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

#ADMIXTURE visualization

cv<-read.csv("admixture_CV.csv")

ggplot(cv)+
  geom_line(aes(K,CV),size=1)+
  #geom_vline(xintercept = 4,linetype="dashed")+
  scale_x_continuous(breaks = c(1:15))+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

tbl1=read.table("P.borealis_stacks.g5mac3dp10ind50g95maf05mdp20p1.recode.1.Q")
tbl1$pop<-pop2
tbl1<-tbl1[c(which(tbl1$pop=="SI18"),which(tbl1$pop=="SU18"),which(tbl1$pop=="SD18"),which(tbl1$pop=="IS21"),which(tbl1$pop=="YS21"),which(tbl1$pop=="UH21"),which(tbl1$pop=="AR18")),]
barplot(t(as.matrix(tbl1)), col=c("blue","red"),xlab=NA, ylab="Ancestry", border=NA,cex.names = 0.75,las=2)

tbl2=read.table("P.borealis_stacks.g5mac3dp10ind50g95maf05mdp20p1.recode.2.Q")
tbl2$pop<-pop2
tbl2$ind<-samp2
row.names(tbl2)<-samp2
tbl2<-tbl2[order(tbl2[,3],tbl2[,1],tbl2[,2]),]
tbl2<-tbl2[c(which(tbl2$pop=="AR18"),which(tbl2$pop=="SI18"),which(tbl2$pop=="SU18"),which(tbl2$pop=="IS21"),which(tbl2$pop=="SD18"),which(tbl2$pop=="YS21"),which(tbl2$pop=="UH21")),]
barplot(t(as.matrix(tbl2)), col=c("red","blue"),xlab=NA, ylab="Ancestry", border=NA,las=2,cex.lab=1.5,cex.names = .75,cex.axis = 1.5)

Qmean<-data.frame(aggregate(V1~pop,tbl2,FUN=mean),V2=aggregate(V2~pop,tbl2,FUN=mean)$V2)
Qmean<-data.frame(V1=Qmean$V1,V2=Qmean$V2,pop=c("AR","S3","S4","S1","S2","OS","S5"))
row.names(Qmean)<-Qmean$pop
Qmean<-Qmean[c(which(Qmean$pop=="AR"),which(Qmean$pop=="S1"),which(Qmean$pop=="S2"),which(Qmean$pop=="S3"),which(Qmean$pop=="S4"),which(Qmean$pop=="S5"),which(Qmean$pop=="OS")),]
Qmean$pop<-factor(Qmean$pop,level=Qmean$pop)

ggplot(data = Qmean)+
  geom_col(aes(x=pop,y=V2),fill="orange")+
  #geom_line(inherit.aes = F,aes(x=Qmean$pop,y=c(NA,Qmean$V2[2],Qmean$V2[3],Qmean$V2[4],Qmean$V2[5],Qmean$V2[6],Qmean$V2[7])),group=1) +
  geom_line(inherit.aes = F,aes(x=Qmean$pop,y=c(NA,NA,NA,Qmean$V2[4],sum(Qmean$V2[4],(Qmean$V2[6]-Qmean$V2[4])/2),Qmean$V2[6],Qmean$V2[7])),group=1,size=1.5) +
  xlab("")+
  ylab("Mean outshore ancestry assignment")+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))


indmiss<-read.table("postfiltmiss.imiss",sep="\t",header = T)

indmiss$Sample<-samp2
indmiss$Pop<-pop2
indmiss<-indmiss[order(indmiss$F_MISS,decreasing = T),]
indmiss$Sample<-factor(indmiss$Sample,levels = indmiss$Sample)
ggplot(data=indmiss)+
  geom_col(aes(x = Sample, y = (F_MISS)))+
  ylab("Fraction of missing SNPs")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.1,hjust = 0.1))

aggregate(N_MISS~Pop,indmiss,FUN=max)
aggregate(N_MISS~Pop,indmiss,FUN=min)
aggregate(F_MISS~Pop,indmiss,FUN=mean)

cv18<-read.csv("admixture_CV_2018.csv")

ggplot(cv18)+
  geom_line(aes(K,CV),size=1)+
  #geom_vline(xintercept = 4,linetype="dashed")+
  scale_x_continuous(breaks = c(1:15))+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

cv21<-read.csv("admixture_CV_2021.csv")

ggplot(cv21)+
  geom_line(aes(K,CV),size=1)+
  #geom_vline(xintercept = 4,linetype="dashed")+
  scale_x_continuous(breaks = c(1:15))+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

tbl2_18=read.table("P.borealis_2018.2.Q")
tbl2_18$pop<-pop2[1:36]
tbl2_18$ind<-samp2[1:36]
row.names(tbl2_18)<-samp2[1:36]
tbl2_18<-tbl2_18[order(tbl2_18[,3],tbl2_18[,1],tbl2_18[,2]),]
tbl2_18<-tbl2_18[c(which(tbl2_18$pop=="SI18"),which(tbl2_18$pop=="SU18"),which(tbl2_18$pop=="SD18"),which(tbl2_18$pop=="AR18")),]
par(mar=c(8,4.4,1,0.5))
barplot(t(as.matrix(tbl2_18)), col=c("blue","red"),xlab=NA, ylab="Ancestry", border=NA,las=2,cex.lab=1.5,cex.names = 1.5,cex.axis = 1.5)

tbl2_21=read.table("P.borealis_2021.2.Q")
tbl2_21$pop<-pop2[37:68]
tbl2_21$ind<-samp2[37:68]
row.names(tbl2_21)<-samp2[37:68]
tbl2_21<-tbl2_21[order(tbl2_21[,3],tbl2_21[,1],tbl2_21[,2]),]
tbl2_21<-tbl2_21[c(which(tbl2_21$pop=="IS21"),which(tbl2_21$pop=="YS21"),which(tbl2_21$pop=="UH21")),]
par(mar=c(8,4.4,1,0.5))
barplot(t(as.matrix(tbl2_21)), col=c("blue","red"),xlab=NA, ylab="Ancestry", border=NA,las=2,cex.lab=1.5,cex.names = 1.5,cex.axis = 1.5)


#Fjarlægð
library(geo)
pos1 <- list(lat = c(65, 66), lon = c(-19, -20))
pos2 <- list(lat = c(64, 65), lon = c(-19, -20))

dists <- arcdist(pos1, pos2)         # pos1 and pos2 are lists of coordinates

std<-read.csv("stodvar.csv")

AR18 <- list(lat = c(as.numeric(std$kastad_breidd)[2], as.numeric(std$hift_breidd)[2]), lon = c(as.numeric(std$kastad_lengd)[2], as.numeric(std$hift_lengd)[2]))
SI18<- list(lat = c(as.numeric(std$kastad_breidd)[4], as.numeric(std$hift_breidd)[4]), lon = c(as.numeric(std$kastad_lengd)[4], as.numeric(std$hift_lengd)[4]))
SU18<- list(lat = c(as.numeric(std$kastad_breidd)[3], as.numeric(std$hift_breidd)[3]), lon = c(as.numeric(std$kastad_lengd)[3], as.numeric(std$hift_lengd)[3]))
IS21<- list(lat = c(as.numeric(std$kastad_breidd)[7], as.numeric(std$hift_breidd)[7]), lon = c(as.numeric(std$kastad_lengd)[7], as.numeric(std$hift_lengd)[7]))
SD18<- list(lat = c(as.numeric(std$kastad_breidd)[1], as.numeric(std$hift_breidd)[1]), lon = c(as.numeric(std$kastad_lengd)[1], as.numeric(std$hift_lengd)[1]))
YS21<- list(lat = c(as.numeric(std$kastad_breidd)[6], as.numeric(std$hift_breidd)[6]), lon = c(as.numeric(std$kastad_lengd)[6], as.numeric(std$hift_lengd)[6]))
UH21 <- list(lat = c(as.numeric(std$kastad_breidd)[5], as.numeric(std$hift_breidd)[5]), lon = c(as.numeric(std$kastad_lengd)[5], as.numeric(std$hift_lengd)[5]))

DIST<-matrix(data=c(
  0,
  mean(arcdist(SI18,SU18)),
  mean(arcdist(SI18,SD18)),
  mean(arcdist(SI18,IS21)),
  mean(arcdist(SI18,YS21)),
  mean(arcdist(SI18,UH21)),
  mean(arcdist(SI18,AR18)),
  #
  mean(arcdist(SU18,SI18)),
  0,
  mean(arcdist(SU18,SD18)),
  mean(arcdist(SU18,IS21)),
  mean(arcdist(SU18,YS21)),
  mean(arcdist(SU18,UH21)),
  mean(arcdist(SU18,AR18)),
  #
  mean(arcdist(SD18,SI18)),
  mean(arcdist(SD18,SU18)),
  0,
  mean(arcdist(SD18,IS21)),
  mean(arcdist(SD18,YS21)),
  mean(arcdist(SD18,UH21)),
  mean(arcdist(SD18,AR18)),
  #
  mean(arcdist(IS21,SI18)),
  mean(arcdist(IS21,SU18)),
  mean(arcdist(IS21,SD18)),
  0,
  mean(arcdist(IS21,YS21)),
  mean(arcdist(IS21,UH21)),
  mean(arcdist(IS21,AR18)),
  #
  mean(arcdist(YS21,SI18)),
  mean(arcdist(YS21,SU18)),
  mean(arcdist(YS21,SD18)),
  mean(arcdist(YS21,IS21)),
  0,
  mean(arcdist(YS21,UH21)),
  mean(arcdist(YS21,AR18)),
  #
  mean(arcdist(UH21,SI18)),
  mean(arcdist(UH21,SU18)),
  mean(arcdist(UH21,SD18)),
  mean(arcdist(UH21,IS21)),
  mean(arcdist(UH21,YS21)),
  0,
  mean(arcdist(UH21,AR18)),
  #
  mean(arcdist(AR18,SI18)),
  mean(arcdist(AR18,SU18)),
  mean(arcdist(AR18,SD18)),
  mean(arcdist(AR18,IS21)),
  mean(arcdist(AR18,YS21)),
  mean(arcdist(AR18,UH21)),
  0),nrow=7)
  
colnames(DIST)<-c("SI18","SU18","SD18","IS21","YS21","UH21","AR18")
  
FST<-read.csv("FST_WC_weighted.csv")
FST<-as.matrix(FST[,-1])

#Öll gögn
vegan::mantel(FST,DIST,permutations=999)
ade4::mantel.rtest(as.dist(FST),as.dist(DIST),)
vegan::mantel(DIST,FST,permutations=999)
ade4::mantel.rtest(as.dist(DIST),as.dist(FST))

#2018 - án Arnarfjarðar
vegan::mantel(FST[c(1,2,3),c(1,2,3)],DIST[c(1,2,3),c(1,2,3)],permutations=999)
vegan::mantel(DIST[c(1,2,3),c(1,2,3)],FST[c(1,2,3),c(1,2,3)],permutations=999)

#2021
vegan::mantel(FST[c(4,5,6),c(4,5,6)],DIST[c(4,5,6),c(4,5,6)],permutations=999)
vegan::mantel(DIST[c(4,5,6),c(4,5,6)],FST[c(4,5,6),c(4,5,6)],permutations=999)

#2018 - öll
vegan::mantel(FST[c(1,2,3,7),c(1,2,3,7)],DIST[c(1,2,3,7),c(1,2,3,7)],permutations=999)

#Öll gögn - á Arnarfjarðar
vegan::mantel(FST[-7,-7],DIST[-7,-7],method="kendall",permutations=999)
vegan::mantel(DIST[-7,-7],FST[-7,-7],method="kendall",permutations=999)
ade4::mantel.rtest(as.dist(FST[-7,-7]),as.dist(DIST[-7,-7]))
ade4::mantel.rtest(as.dist(DIST[-7,-7]),as.dist(FST[-7,-7]))

plot(as.dist(DIST),as.dist(FST))
plot(as.dist(DIST[c(1,2,3),c(1,2,3)]),as.dist(FST[c(1,2,3),c(1,2,3)]))
plot(as.dist(DIST[c(1,2,3,7),c(1,2,3,7)],as.dist(FST[c(1,2,3,7),c(1,2,3,7)])))
plot(as.dist(DIST[c(4,5,6),c(4,5,6)]),as.dist(FST[c(4,5,6),c(4,5,6)]))

#Fjarlægð sem sýnir ekki IBD 
plot(log(as.dist(DIST[-7,-7])),as.dist(FST[-7,-7]/(1-FST[-7,-7])))

arcdist(lat=AR18[[1]][1],lon=AR18[[2]][1],lat1=AR18[[1]][2],lon1=AR18[[2]][2])

env<-std[,c(6,7,8,9,12,13,14,16)]
row.names(env)<-c("SU18","AR18","SI18","SD21","YS21","IS21","UH21")
DIST
ad<-vegan::adonis2(as.dist(FST)~kastad_breidd*kastad_lengd*hift_breidd*hift_lengd*botnhiti*yfirbordshiti*lofthiti*toglengd,data=env, permutations = 99)
summary(ad)

require(hierfstat)
?hierfstat
#prepare sample for OutFLANK analysis
G2 <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno),dimnames = list(pos_loc,colnames(geno)))

G2[geno %in% c("0/0")] <- 00
G2[geno %in% c("0/1", "1/0")] <- 10
G2[geno %in% c("1/1")] <- 11

dat<-data.frame(t(G2))
dat$Pop<-pop2
dat<-dat[,c(1472,1:1471)]
betas(dat,nboot=1000,lim=c(0.025,0.975),diploid=TRUE,betaijT=FALSE)
boot.ppbetas(dat=dat,nboot=100,quant=c(0.025,0.975),diploid=TRUE,digits=4)

require(dartR)
