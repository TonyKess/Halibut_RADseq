#rm(list = ls())
#Quick script for inversion. (modified from S. Lehnert)
#Need bed/bim/fam file from plink with chromosome position information included
#Can have whole genome and indicate which chromosome to run
#Seems to work better with all data (not filtered for maf 0.05)

###Inversion packages
library(snpStats)
library(inveRsion)
install.packages("inveRsion")
#Package for changing characters
library(stringr) 

#Working directory
setwd("/Users/bradburygenomicslab/Desktop/Projects/Halibut_FinalRAD/")
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV  --chr 15 --out Halibut_CHR15 --recode --make-bed")

#Note chromosome names have to be NUMBERS* 
data1 <- snpStats::read.plink(bed = "Halibut_CHR15.bed",
                              bim = "Halibut_CHR15.bim",
                              fam = "Halibut_CHR15.fam")

#Set up genotype data for Chr 15 
ThingCHR15 <- setUpGenoDatSNPmat(Chr = "15", Geno = data1$genotypes, Annot = data1$map)

#Run hapcode  - Can change minAllele (this is minimum frequecny for the inversion)
hapCode_thing15 <- codeHaplo(ThingCHR15, minAllele = 0.1, saveRes = TRUE) #Default block size of 3 SNPs used
window<-1 #1 MBp - can change

#Scan for inversion from hapCode data
scanRes_thing15 <- scanInv(hapCode_thing15 , window = window, geno=TRUE, saveRes=TRUE)

#Get list of inversions with breakpoints
invList_thing15<- listInv(scanRes_thing15 , hapCode = hapCode_thing15, geno=TRUE, all = F, saveRes = TRUE, thBic=0, saveROI = TRUE)

#Gives identified inversions -- look for one with highest number of models and high MaxBIC
invList_thing15

#These steps above take a while to run so I would save results if need to load later..
save.image("CHR15_inversiondata.RData")

##Choose inversion that you want to look at (wROI = row # in invList that looks to have the best support - region of interest)
#Get assignments to karyotype class
r15 <- getClassif(invList_thing15, geno=TRUE, wROI = 1)

#Merge classification with your individual info (from fam file)
individuals=read.table("Halibut_CHR15.fam",header=F)
individuals

inv_assignments=cbind(r15, individuals)
length(inv_assignments$hom[inv_assignments$hom==1])
length(inv_assignments$hom[inv_assignments$homInv==1])
length(inv_assignments$het[inv_assignments$het==1])
