library(pcadapt)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggman)
library(qvalue)
setwd("~/Desktop/Projects/Halibut_FinalRAD/")

#Check basic sample info####
FAM <- read.delim("Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV.fam", stringsAsFactors = F, header = F, sep = "")
Metadata <- read.delim("Metadata_Halibut_ENV_noNAnoWF.txt", header = T, stringsAsFactors = F)
ID <- FAM[1]
colnames(ID)[1] <- "ID"
Metadata_sorted <- inner_join(ID, Metadata)


#pcadapt on entire dataset####
Halib<- read.pcadapt("Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV.ped", type = "ped")
PCs <- pcadapt(Halib, K = 3, min.maf = 0.01)
Metadata_sorted <- inner_join(ID, Metadata)
PC_scores_pop <- as.data.frame(cbind(Metadata_sorted, PCs$scores))
cor.test(PC_scores_pop$'3', PC_scores_pop$TEMP_PC)

plot(PCs, option = "manhattan" ) + theme_classic()
plot(PCs, option = "screeplot" ) + theme_classic()
#plot multiple axes
plot(PCs, option="scores", pop =  Metadata_sorted$Sex, i = 1, j =2 ) + theme_classic() 
plot(PCs, option="scores", pop =  Metadata_sorted$Region, i = 1, j = 3) + theme_classic()

#Plot loadings 
Chrom_map <- read.delim("Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV.map", stringsAsFactors = F, header = F)
colnames(Chrom_map) <- c("Chrom", "SNP", "CM", "BP")
PVALS <- PCs$pvalues
PCMAP <- as.data.frame(cbind(Chrom_map, PVALS))
QVALS <- qvalue(PCMAP$PVALS)
PCMAP$QVALS <- QVALS$qvalues
#Get outliers
OL_K3_Sex_Included <- PCMAP[which(PCMAP$QVALS<0.05),]
PCA_Man <- ggman(PCMAP, chrom = "Chrom", pvalue = "QVALS", snp = "SNP", bp="BP", pointSize = 1, title = "Halibut", xlabel = "Chromosome", ylabel = "-log10(q value)" ,sigLine = -log10(0.05)  ) + theme_classic()
ggmanHighlight(PCA_Man, OL_K3_Sex_Included$SNP)
OL_K3_Autosome <- OL_K3_Sex_Included[which(OL_K3_Sex_Included$Chrom!="18"),]
OL_K3_NoInversions_sex <- OL_K3_Autosome[which(OL_K3_Autosome$Chrom!="15"),] 

write.table(OL_K3_NoInversions_sex, "OL_K3_NoInversions_sex.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(OL_K3_Autosome, "OL_K3_Autosome_Outliers.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#subsample to prevent illustrator from exploding
Subsamp <-sample_n(Chrom_map, 15000 )
Subsamp_PC <- inner_join(Subsamp, PCMAP)
Subsamp_PC_withOL <- distinct(rbind(Subsamp_PC, OL_K3_Sex_Included))

#plot by qvalue
PCA_Man <- ggman(Subsamp_PC_withOL, chrom = "Chrom", pvalue = "QVALS", snp = "SNP", bp="BP", pointSize = 1, title = "Halibut", xlabel = "Chromosome", ylabel = "-log10(q value)" ,sigLine = -log10(0.05)  ) + theme_classic()
ggmanHighlight(PCA_Man, OL_K3_Sex_Included$SNP)

#per-axis loadings...####
PCs_compwise <- pcadapt(Halib, K = 3, min.maf = 0.01, method = "componentwise")
PC_scores_compwise_pop <- as.data.frame(cbind(Metadata_sorted, PCs_compwise$scores))

PVALS_PCs_compwise <- PCs_compwise$pvalues
PCMAP_compwise <- as.data.frame(cbind(Chrom_map, PVALS_PCs_compwise))
colnames(PCMAP_compwise)[5:7] <- c("PC1", "PC2", "PC3")

#get OL
QVALS_PC1 <- qvalue(PCMAP_compwise$PC1)
PCMAP_compwise$PC1_QVALS <- QVALS_PC1$qvalues
OL_PC1 <- PCMAP_compwise[which(PCMAP_compwise$PC1_QVALS <0.05),]


QVALS_PC2 <- qvalue(PCMAP_compwise$PC2)
PCMAP_compwise$PC2_QVALS <- QVALS_PC2$qvalues
OL_PC2 <- PCMAP_compwise[which(PCMAP_compwise$PC2_QVALS <0.05),]

QVALS_PC3 <- qvalue(PCMAP_compwise$PC3)
PCMAP_compwise$PC3_QVALS <- QVALS_PC3$qvalues
OL_PC3 <- PCMAP_compwise[which(PCMAP_compwise$PC3_QVALS <0.05),]


Subsamp_PC_compwise <- inner_join(Subsamp, PCMAP_compwise)
Subsamp_PC_compwise_withOL <- distinct(rbind(Subsamp_PC_compwise, OL_PC1, OL_PC2, OL_PC3))

OL_compwise <- distinct(rbind(OL_PC1, OL_PC2, OL_PC3))
inner_join(OL_compwise, OL_K3_Sex_Included, by = "SNP")
anti_join(OL_compwise, OL_K3_Sex_Included, by = "SNP")

length(OL_compwise$SNP[OL_compwise$SNP %in% OL_K3_Sex_Included$SNP])/length(unique(c(OL_compwise$SNP, OL_K3_Sex_Included$SNP)))

#plot
PCA_Man_PC1 <- ggman(Subsamp_PC_compwise_withOL, chrom = "Chrom", pvalue = "PC1_QVALS", snp = "SNP", bp="BP", pointSize = 1, title = "Halibut", xlabel = "Chromosome", ylabel = "-log10(q value)" ,sigLine = -log10(0.05), ymax = 10) + theme_classic()
ggmanHighlight(PCA_Man_PC1, OL_PC1$SNP)

PCA_Man_PC2 <- ggman(Subsamp_PC_compwise_withOL, chrom = "Chrom", pvalue = "PC2_QVALS", snp = "SNP", bp="BP", pointSize = 1, title = "Halibut", xlabel = "Chromosome", ylabel = "-log10(q value)" ,sigLine = -log10(0.05), ymax = 10  ) + theme_classic()
ggmanHighlight(PCA_Man_PC2, OL_PC2$SNP)

PCA_Man_PC3 <- ggman(Subsamp_PC_compwise_withOL, chrom = "Chrom", pvalue = "PC3_QVALS", snp = "SNP", bp="BP", pointSize = 1, title = "Halibut", xlabel = "Chromosome", ylabel = "-log10(q value)" ,sigLine = -log10(0.05), ymax = 10) + theme_classic()
ggmanHighlight(PCA_Man_PC3, OL_PC3$SNP)
OL_PC3$SNP   %in% Subsamp_PC_compwise_withOL


#remove CHR15 and CHR18 and do PCA again
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV --out Halibut_genome2_d15gt90HWE_NOCHR1518 --not-chr 18, 15 --recode --make-bed --chr-set 24 ")
Halib_noinv_nosex<- read.pcadapt("Halibut_genome2_d15gt90HWE_NOCHR1518.ped", type = "ped")

PCs_NINS_compwise <- pcadapt(Halib_noinv_nosex, K = 3, min.maf = 0.01, method = "componentwise")
PC_NINS_scores_compwise_pop <- as.data.frame(cbind(Metadata_sorted, PCs_NINS_compwise$scores))

#plot multiple axes
plot(PCs_NINS_compwise, option="scores", pop =  Metadata_sorted$Region, i = 1, j =2 ) + theme_classic() 
#Get inversion haplotypes
#Plot loadings 
Chrom_map_NINS <- read.delim("Halibut_genome2_d15gt90HWE12_NOCHR1518.map", stringsAsFactors = F, header = F)
colnames(Chrom_map_NINS) <- c("Chrom", "SNP", "CM", "BP")

PVALS_PC_NINS_compwise <- PCs_NINS_compwise$pvalues
PCMAP_NINS_compwise <- as.data.frame(cbind(Chrom_map_NINS, PVALS_PC_NINS_compwise))
colnames(PCMAP_NINS_compwise)[5:7] <- c("PC1", "PC2", "PC3")

#get OL
QVALS_NINS_PC1 <- qvalue(PCMAP_NINS_compwise$PC1)
PCMAP_NINS_compwise$PC1_QVALS <- QVALS_NINS_PC1$qvalues
OL_NINS_PC1 <- PCMAP_NINS_compwise[which(PCMAP_NINS_compwise$PC1_QVALS <0.05),]


Subsamp_NINS_PC_compwise <- inner_join(Subsamp[1:4], PCMAP_NINS_compwise, )
Subsamp_NINS_PC_compwise_withOL <- distinct(rbind(Subsamp_NINS_PC_compwise, OL_NINS_PC1))

PCA_Man_PC1 <- ggman(Subsamp_NINS_PC_compwise_withOL, chrom = "Chrom", pvalue = "PC1_QVALS", snp = "SNP", bp="BP", pointSize = 1, title = "Halibut", xlabel = "Chromosome", ylabel = "-log10(q value)" ,sigLine = -log10(0.05), ymax = 10) + theme_classic()
ggmanHighlight(PCA_Man_PC1, OL_NINS_PC1$SNP)
OL_NINS_PC1$SNP   %in% Subsamp_NINS_PC_compwise_withOL$SNP

OL_compwise <- unique(c(OL_PC1$SNP, OL_PC2$SNP, OL_NINS_PC1$SNP))

#Analyze PCA/SNMF with LD pruning ####
#LD PRUNE 
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_genome2_d15gt90HWE_NOCHR1518 --indep-pairwise 50 5 0.5 --out Halibut_window_LD --chr-set 24 ")
#Plink filter
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_genome2_d15gt90HWE12_NOCHR1518 --exclude Halibut_window_LD.prune.out --out Halibut_genome2_d15gt90HWE12_NOCHR1518_LE --recode --make-bed --chr-set 24 ")
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_genome2_d15gt90HWE12_NOCHR1518 --exclude Halibut_window_LD.prune.out --out Halibut_genome2_d15gt90HWE12_NOCHR1518_LE12  --recode12 --make-bed --chr-set 24 ")
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_genome2_d15gt90HWE_NOCHR1518  --out Halibut_genome2_d15gt90HWE_NOCHR1518_LE --recode vcf --make-bed --chr-set 24 ")

#Get RAD locus info from SNP ids in .map file (update with sed)
RADloc_strand <- fread("Halibut_Chrom_POS_Locus_LP_strand.txt", data.table = F, stringsAsFactors = F)
colnames(RADloc_strand) <- c("Chrom", "BP", "RADLoc", "RADLoc_Pos", "Strand")
Chrom_RAD_map <- inner_join(Chrom_map, RADloc_strand)

Chrom_map_LE <- fread("Halibut_genome2_d15gt90HWE12_NOCHR1518_LE.map", data.table = F, stringsAsFactors = F)
colnames(Chrom_map_LE) <- c("Chrom", "SNP", "CM", "BP")
Chrom_RAD_map_LE <- inner_join(Chrom_map_LE, RADloc_strand)

#pcadapt on LE
Halib_LE<- read.pcadapt("Halibut_genome2_d15gt90HWE12_NOCHR1518_LE.ped", type = "ped")
PCs_LE <- pcadapt(Halib_LE, K = 3, min.maf = 0.01)
PCs_LE$singular.values ^2 * 100
Metadata_sorted <- inner_join(ID, Metadata)
PC_LE_scores_pop <- as.data.frame(cbind(Metadata_sorted, PCs_LE$scores))

plot(PCs_LE, option = "manhattan" ) + theme_classic()
plot(PCs_LE, option = "screeplot" ) + theme_classic()
#plot multiple axes
plot(PCs_LE, option="scores", pop =  Metadata_sorted$Region, i = 1, j =2 ) + theme_classic() 
#Get inversion haplotypes
#Plot loadings 

#Single RAD ####
Single_RAD <- Chrom_RAD_map$SNP[!duplicated(Chrom_RAD_map$RADLoc)]
#remove one variant that somehow escaped duplicate
Single_RAD <- Single_RAD[!(Single_RAD %in% "378936:81:+")]

write.table(Single_RAD, "Single_RAD.txt", sep = "\t",col.names = F, row.names = F, quote = F)

system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_genome2_d15gt90HWE_NOCHR1518 --extract Single_RAD.txt --out Halibut_genome2_d15gt90HWE12_NOCHR1518_singleRAD --recode --make-bed --chr-set 24 ")
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_genome2_d15gt90HWE_NOCHR1518 --extract Single_RAD.txt --out Halibut_genome2_d15gt90HWE12_NOCHR1518_singleRAD12 --recode12 --make-bed --chr-set 24 ")
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_genome2_d15gt90HWE_NOCHR1518 --extract Single_RAD.txt --out Halibut_genome2_d15gt90HWE12_NOCHR1518_singleRAD --recode vcf --make-bed --chr-set 24 ")


Halib_SR<- read.pcadapt("Halibut_genome2_d15gt90HWE12_NOCHR1518_singleRAD.ped", type = "ped")
PCs_SR <- pcadapt(Halib_SR, K = 3, min.maf = 0.01)
PC_SR_scores_pop <- as.data.frame(cbind(Metadata_sorted, PCs_SR$scores))
(PCs_SR$singular.values ^ 2) * 100
#plot multiple axes
plot(PCs_SR, option="scores", pop =  Metadata_sorted$Region, i = 1, j =2 ) + theme_classic() 


#pcadapt on Atlantic ####
system("~/Desktop/Software/plink_mac_20200219/plink --file  Halibut_genome2_d15gt90HWE12_NOCHR1518_LE  --keep Notgulf_inds_plink.txt --out Halibut_genome2_d15gt90HWE12_NOCHR1518_notGulf --recode --make-bed --chr-set 24 ")
Halib_Atlantic<- read.pcadapt("Halibut_genome2_d15gt90HWE12_NOCHR1518_notGulf.ped", type = "ped")
PCs_Atl <- pcadapt(Halib_Atlantic, K = 3, min.maf = 0.01)

Atl_FAM <- read.delim("Halibut_genome2_d15gt90HWE12_NOCHR1518_notGulf.fam", stringsAsFactors = F, header = F, sep = "")
Atl_ID <- Atl_FAM[1]
colnames(Atl_ID)[1] <- "ID"
Atl_Metadata_sorted <- inner_join(Atl_ID, Metadata)

PC_Atl_scores_pop <- as.data.frame(cbind(Atl_Metadata_sorted, PCs_Atl$scores))

plot(PCs_Atl, option = "manhattan" ) + theme_classic()
plot(PCs_Atl, option = "screeplot" ) + theme_classic()
#plot multiple axes
plot(PCs_Atl, option="scores", pop =Atl_Metadata_sorted$Region, i = 1, j =2 ) + theme_classic() 

#pcadapt on Gulf ####

system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_genome2_d15gt90HWE12_NOCHR1518_LE -keep Gulf_inds_plink.txt --out Halibut_genome2_d15gt90HWE12_NOCHR1518_Gulf --recode --make-bed --chr-set 24 ")
Halib_Gulf<- read.pcadapt("Halibut_genome2_d15gt90HWE12_NOCHR1518_Gulf.ped", type = "ped")
PCs_Gulf <- pcadapt(Halib_Gulf, K = 3, min.maf = 0.01)

Gulf_FAM <- read.delim("Halibut_genome2_d15gt90HWE12_NOCHR1518_Gulf.fam", stringsAsFactors = F, header = F, sep = "")
Gulf_ID <- Gulf_FAM[1]
colnames(Gulf_ID)[1] <- "ID"
Gulf_Metadata_sorted <- inner_join(Gulf_ID, Metadata)

PC_Gulf_scores_pop <- as.data.frame(cbind(Gulf_Metadata_sorted, PCs_Gulf$scores))

plot(PCs_Atl, option = "manhattan" ) + theme_classic()
plot(PCs_Atl, option = "screeplot" ) + theme_classic()
#plot multiple axes
plot(PCs_Gulf, option="scores", pop =Gulf_Metadata_sorted$Region, i = 1, j =2 ) + theme_classic() 
#Get inversion haplotypes
#Plot loadings 

library(LEA)
library(lfmm)
###SNMF #ALL  ####
ped2lfmm(input.file = "Halibut_genome2_d15gt90HWE12_NOCHR1518.ped")
pc = pca("Halibut_genome2_d15gt90HWE12_NOCHR1518.lfmm") 
tc = tracy.widom(pc)
plot(tc$percentage)
project2 = NULL
project2 = snmf(input.file = "Halibut_genome2_d15gt90HWE12_NOCHR1518.lfmm", K = 1:5, entropy = TRUE, project = "new", CPU = 16 )
plot(project2, col = "blue", pch = 19, cex = 1.2)
plot(project2)
best = which.min(cross.entropy(project2, K = 2))
my.colors <- c( "dodgerblue", "red", "purple")
bg2 <- colorRampPalette(c("dodgerblue", "red", "purple", "cyan", "orange", "goldenrod", "orchid"))(n=35)
barchart(project2, K = 2, run = best,
         border = NA, space = 0,
         col = bg2,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = ID_REGION$Region[bp$order], las=1,
     cex.axis = .3)

library(tidyverse)
POP <- data.frame(cbind(Metadata$ID, Metadata$Region), stringsAsFactors = F)
colnames(POP) <- c("ID", "Regions")
HBUT_IDs <-as.data.frame(system("awk '{print $2}' Halibut_genome2_d15gt90HWE12_NOCHR1518.ped", intern = T), stringsAsFactors = F)

colnames(HBUT_IDs) <- "ID"
ID_REGION <- dplyr::inner_join(HBUT_IDs, POP)
ID_REGION$ID

Qval <- data.frame(cbind(ID_REGION, Q(object = project2, K = 2, run = 1)), stringsAsFactors = F)
Qval_ord <- Qval[order(Qval$Region),]
tbl = Qval_ord
rownames(tbl) <- Qval_ord$ID

plot_data <-  tbl %>% 
  gather('pop', 'prob', V1:V2) %>% 
  group_by(ID)

ggplot(plot_data, aes(ID, prob, fill = pop)) +
  geom_col() +
  facet_grid(~Regions, scales = 'free', space = 'free')

wilcox.test(Qval_ord$V1[Qval_ord$Regions %in% "Gulf of St. Lawrence"], Qval_ord$V1[!(Qval_ord$Regions %in% "Gulf of St. Lawrence")])


###SNMF #Linkage equilib ####
ped2lfmm(input.file = "Halibut_genome2_d15gt90HWE12_NOCHR1518_LE12.ped")
pc = pca("Halibut_genome2_d15gt90HWE12_NOCHR1518_LE12.lfmm") 
tc = tracy.widom(pc)
plot(tc$percentage)
project2 = NULL
project2 = snmf(input.file = "Halibut_genome2_d15gt90HWE12_NOCHR1518_LE12.lfmm", K = 1:5, entropy = TRUE, project = "new", CPU = 16 )
plot(project2, col = "blue", pch = 19, cex = 1.2)

best = which.min(cross.entropy(project2, K = 2))
my.colors <- c( "dodgerblue", "red", "purple")
bg2 <- colorRampPalette(c("dodgerblue", "red", "purple", "cyan", "orange", "goldenrod", "orchid"))(n=35)
barchart(project2, K = 2, run = best,
         border = NA, space = 0,
         col = bg2,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = ID_REGION$Region[bp$order], las=1,
     cex.axis = .3)

POP <- data.frame(cbind(Metadata$ID, Metadata$Region), stringsAsFactors = F)
colnames(POP) <- c("ID", "Regions")
HBUT_IDs <-as.data.frame(system("awk '{print $2}' Halibut_genome2_d15gt90HWE12_NOCHR1518.ped", intern = T), stringsAsFactors = F)

colnames(HBUT_IDs) <- "ID"
ID_REGION <- dplyr::inner_join(HBUT_IDs, POP)
ID_REGION$ID

Qval <- data.frame(cbind(ID_REGION, Q(object = project2, K = 2, run = 1)), stringsAsFactors = F)
Qval_ord <- Qval[order(Qval$Region),]
tbl = Qval_ord
rownames(tbl) <- Qval_ord$ID

plot_data <-  tbl %>% 
  gather('pop', 'prob', V1:V2) %>% 
  group_by(ID)

ggplot(plot_data, aes(ID, prob, fill = pop)) +
  geom_col() +
  facet_grid(~Regions, scales = 'free', space = 'free')

wilcox.test(Qval_ord$V1[Qval_ord$Regions %in% "Gulf of St. Lawrence"], Qval_ord$V1[!(Qval_ord$Regions %in% "Gulf of St. Lawrence")])


#single rad SMNF ####
ped2lfmm(input.file = "Halibut_genome2_d15gt90HWE12_NOCHR1518_singleRAD12.ped")
pc = pca("Halibut_genome2_d15gt90HWE12_NOCHR1518_singleRAD12.lfmm") 
tc = tracy.widom(pc)
plot(tc$percentage)
project2 = NULL
project2 = snmf(input.file = "Halibut_genome2_d15gt90HWE12_NOCHR1518_singleRAD12.lfmm", K = 1:5, entropy = TRUE, project = "new", CPU = 16 )
plot(project2, col = "blue", pch = 19, cex = 1.2)

best = which.min(cross.entropy(project2, K = 2))
my.colors <- c( "dodgerblue", "red", "purple")
bg2 <- colorRampPalette(c("dodgerblue", "red", "purple", "cyan", "orange", "goldenrod", "orchid"))(n=35)
barchart(project2, K = 2, run = best,
         border = NA, space = 0,
         col = bg2,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = ID_REGION$Region[bp$order], las=1,
     cex.axis = .3)


Qval <- data.frame(cbind(ID_REGION, Q(object = project2, K = 2, run = 1)), stringsAsFactors = F)
Qval_ord <- Qval[order(Qval$Region),]
tbl = Qval_ord
rownames(tbl) <- Qval_ord$ID

plot_data <-  tbl %>% 
  gather('pop', 'prob', V1:V2) %>% 
  group_by(ID)

ggplot(plot_data, aes(ID, prob, fill = pop)) +
  geom_col() +
  facet_grid(~Regions, scales = 'free', space = 'free')

wilcox.test(Qval_ord$V1[Qval_ord$Regions %in% "Gulf of St. Lawrence"], Qval_ord$V1[!(Qval_ord$Regions %in% "Gulf of St. Lawrence")])
