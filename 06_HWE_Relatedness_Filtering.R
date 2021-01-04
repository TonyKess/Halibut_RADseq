#Use pcadapt to explore clustering etc.

library(pcadapt)
Halib<- read.pcadapt("Halibut_finalgenome_mindp15_24chrom_ming90.ped", type = "ped")
PCs <- pcadapt(Halib, K = 2, min.maf = 0.01)
FAM <- read.delim("Halibut_finalgenome_mindp15_24chrom_ming90.fam", stringsAsFactors = F, header = F, sep = "")
Pop <- read.delim("Halibut_metadata_main.txt", header = T, stringsAsFactors = F)
ID <- FAM[1]
colnames(ID)[1] <- "ID"
library(dplyr)
Pop_info <- inner_join(ID, Pop)
library(ggplot2)
PC_scores_pop <- as.data.frame(cbind(Pop_info, PCs$scores))

plot(PCs, option = "manhattan" ) + theme_classic()
plot(PCs, option = "screeplot" ) + theme_classic()
plot(PCs, option="scores", pop = Pop_info$Region) + theme_classic()

#plot additional axes
plot(PCs, option="scores", pop = Pop_info$Region, i = 1, j =2 ) + theme_classic() + geom_hline(yintercept = -0.05) + geom_hline(yintercept = 0.05)

#Plot loadings
Chrom_map <- read.delim("Halibut_finalgenome_mindp15_24chrom_ming90.map", stringsAsFactors = F, header = F)
colnames(Chrom_map) <- c("Chrom", "SNP", "CM", "BP")
PVALS <- PCs$pvalues
PCMAP <- as.data.frame(cbind(Chrom_map, PVALS))
library(ggman)
ggman(PCMAP, chrom = "Chrom", pvalue = "PVALS", snp = "SNP", bp="BP", pointSize = 1, title = "Halibut", xlabel = "Chromosome", ymax = 100 ) + theme_classic()

#remove structured fish
no_structure_fish <- PC_scores_pop$ID[which(PC_scores_pop$`2`> (-0.05) & PC_scores_pop$`2` < (0.05))]
structure_fish1 <- FAM$V1[!(FAM$V1 %in% no_structure_fish)]
structure_fish1 %in% Related_INDS
write.table(cbind(structure_fish1, structure_fish1 ), "structure_fish1.txt", col.names = F, sep = "\t", quote = F, row.names = F )
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90 --remove structure_fish1.txt --out Halibut_finalgenome_mindp15_24chrom_ming90_noabcluster --recode --make-bed --chr-set 24 ")
#do again without those structured fish
library(pcadapt)
Halib<- read.pcadapt("Halibut_finalgenome_mindp15_24chrom_ming90_noabcluster.ped", type = "ped")
PCs <- pcadapt(Halib, K = 5, min.maf = 0.01)
FAM <- read.delim("Halibut_finalgenome_mindp15_24chrom_ming90_noabcluster.fam", stringsAsFactors = F, header = F, sep = "")
Pop <- read.delim("Halibut_metadata_main.txt", header = T, stringsAsFactors = F)
ID <- FAM[1]
colnames(ID)[1] <- "ID"
library(dplyr)
Pop_info <- inner_join(ID, Pop)
library(ggplot2)
PC_scores_pop <- as.data.frame(cbind(Pop_info, PCs$scores))

plot(PCs, option = "manhattan" ) + theme_classic()
plot(PCs, option = "screeplot" ) + theme_classic()
plot(PCs, option="scores", pop = Pop_info$Region) + theme_classic()
#plot additional axes
plot(PCs, option="scores", pop = Pop_info$Region, i = 3, j =4 ) + theme_classic() 
structure_fish2 <- PC_scores_pop$ID[which(PC_scores_pop$`3` > (0.1))]
#And again
write.table(cbind(structure_fish2, structure_fish2), "structure_fish2.txt", col.names = F, sep = "\t", quote = F, row.names = F )

system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_noabcluster --remove structure_fish2.txt --out Halibut_finalgenome_mindp15_24chrom_ming90_noabcluster --recode --make-bed --chr-set 24 ")

#do again with structured fish - 751 samples now
library(pcadapt)
Halib<- read.pcadapt("Halibut_finalgenome_mindp15_24chrom_ming90_noabcluster.ped", type = "ped")
PCs <- pcadapt(Halib, K = 5, min.maf = 0.01)
FAM <- read.delim("Halibut_finalgenome_mindp15_24chrom_ming90_noabcluster.fam", stringsAsFactors = F, header = F, sep = "")
Pop <- read.delim("Halibut_metadata_main.txt", header = T, stringsAsFactors = F)
ID <- FAM[1]
colnames(ID)[1] <- "ID"
library(dplyr)
Pop_info <- inner_join(ID, Pop)
library(ggplot2)
PC_scores_pop <- as.data.frame(cbind(Pop_info, PCs$scores))

plot(PCs, option = "manhattan" ) + theme_classic()
plot(PCs, option = "screeplot" ) + theme_classic()
plot(PCs, option="scores", pop = Pop_info$Region) + theme_classic()
#plot additional axes
plot(PCs, option="scores", pop = Pop_info$Region, i = 1, j =3 ) + theme_classic() 
All_structurefish <- c(structure_fish1, structure_fish2)

#Relatedness_compare
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_noabcluster --genome --out relate_WFout --chr-set 24")
system("cat structure_fish1.txt structure_fish2.txt  > all_structurefish.txt")
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90 --keep  all_structurefish.txt --genome --out relate_WFin --chr-set 24")

Relate_WFout <- fread("relate_WFout.genome", data.table = F, stringsAsFactors = F)
Relate_WFin <- fread("relate_WFin.genome", data.table = F, stringsAsFactors = F)

mean(Relate_WFout$PI_HAT)
mean(Relate_WFin$PI_HAT)
ggplot() + geom_histogram(data = Relate_WFout, aes(x = PI_HAT))
ggplot() + geom_histogram(data = Relate_WFin, aes(x = PI_HAT))

ggplot() + geom_density(data = Relate_WFout, aes(x = PI_HAT, fill = "No abbertant PCA clustering", alpha = 0.5) )+
  geom_density(data = Relate_WFin, aes(x = PI_HAT, fill = "Abbertant PCA clustering", alpha = 0.5))  + theme_classic()

#HWE check
GULF_INDS <- Pop$ID[Pop$Region %in% "Gulf of St. Lawrence"] 
GULF_INDS_FILE <- cbind(GULF_INDS, GULF_INDS)
write.table(GULF_INDS_FILE, "GULF_INDS.txt", sep = "\t", quote = F, col.names = F, row.names = F)

#
#Filter
system("~/Desktop/plink_linux/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_noabcluster  --keep GULF_INDS.txt --hardy --out GULF --chr-set 24 ")
system("~/Desktop/plink_linux/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_noabcluster --remove GULF_INDS.txt --hardy --out NOGULF --chr-set 24 ")


GULF_HWE<- read.csv("GULF.hwe", sep="", stringsAsFactors=FALSE)

#Bonferroni
length(GULF_HWE$SNP[which(GULF_HWE$P<(0.05/length(GULF_HWE$P)))])

#FDR
library(qvalue)
QVALS_GULF <- qvalue(GULF_HWE$P)
GULF_HWE$Q <- QVALS_GULF$qvalues
length(GULF_HWE$SNP[which(GULF_HWE$Q<0.05)])

NOGULF_HWE<- read.csv("NOGULF.hwe", sep="", stringsAsFactors=FALSE)

#Bonferroni
length(NOGULF_HWE$SNP[which(NOGULF_HWE$P<(0.05/length(NOGULF_HWE$P)))])

#FDR
QVALS_NOGULF <- qvalue(NOGULF_HWE$P)
NOGULF_HWE$Q <- QVALS_NOGULF$qvalues
length(NOGULF_HWE$SNP[which(NOGULF_HWE$Q<0.05)])


FDR_HWE_NOGULF <- NOGULF_HWE$SNP[which(NOGULF_HWE$Q<0.05)]
FDR_HWE_GULF <- GULF_HWE$SNP[which(GULF_HWE$Q<0.05)]

#Get out of HWE across comparisons
Shared_OutofHWE <- FDR_HWE_NOGULF[FDR_HWE_NOGULF %in% FDR_HWE_GULF]
Out_of_HWE <- Chrom_map[Chrom_map$SNP %in% Shared_OutofHWE,]

#Exclude those in putative sex chromosome
CHR18_map <- Chrom_map[Chrom_map$Chrom %in% "18",]
Sex_chrom_SNP <- CHR18_map$SNP[which(CHR18_map$BP>1.2e+07)]
Out_of_HWE_Finalized <- Out_of_HWE$SNP[!(Out_of_HWE$SNP %in% Sex_chrom_SNP)]

#Filter
write.table(Out_of_HWE_Finalized, "OutofHWE.txt", col.names = F, sep = "\t", quote = F, row.names = F)
system("~/Desktop/plink_linux/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_noabcluster --exclude OutofHWE.txt --chr-set 24 --recode --make-bed --out Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF ")

#Test what we get back now
Halib<- read.pcadapt("Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF.ped", type = "ped")structured
PCs <- pcadapt(Halib, K = 10, min.maf = 0.01)
FAM <- read.delim("Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF.fam", stringsAsFactors = F, header = F, sep = "")
Metadata <- read.delim("Halibut_metadata_main.txt", header = T, stringsAsFactors = F)
ID <- FAM[1]
colnames(ID)[1] <- "ID"
Metadata_sorted <- inner_join(ID, Metadata)
PC_scores_pop <- as.data.frame(cbind(Metadata_sorted, PCs$scores))

plot(PCs, option = "manhattan" ) + theme_classic()
plot(PCs, option = "screeplot" ) + theme_classic()
plot(PCs, option="scores", pop = Pop_info$Region) + theme_classic()
#plot additional axes
plot(PCs, option="scores", pop = Pop_info$Region, i = 1, j =3 ) + theme_classic() 

#Filter for individuals with complete environmental data too
Metadata_ENV<- read.delim("~/Desktop/Tony/Halib_allgenomealigned/Halibut_Metadata_ENV_Complete.txt", header = T, stringsAsFactors = F)
colnames(ID)
ENV_and_notstructured_and_goodSNPs <- inner_join(ID, Metadata_ENV)
write.table(ENV_and_notstructured_and_goodSNPs, "Metadata_Halibut_ENV_noNAnoWF.txt", col.names = T, sep = "\t", row.names = F, quote = F)
write.table(cbind(ENV_and_notstructured_and_goodSNPs$ID, ENV_and_notstructured_and_goodSNPs$ID), "tokeep_Halibut_ENV_noNAnoWF.txt", col.names = T, sep = "\t", row.names = F, quote = F)

#Make a genotype file with all duplicate/high relatedness individuals removed, no missing geographic info, SNPs with consistent patterns of HWE deviation removed
system("~/Desktop/plink_linux/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF --keep tokeep_Halibut_ENV_noNAnoWF.txt  --out Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV --recode --make-bed --chr-set 24 ")
)
