library(pcadapt)
library(data.table)
library(dplyr)
library(ggplot2)
library(rsed)
library(tidyr)
library(marmap)
library(stringr)
library(maps)
library(lattice)

###PCA#####
setwd("~/Desktop/Projects/Halibut_FinalRAD/")

#Import metadata
FAM <- read.delim("Halibut_CHR15.fam", stringsAsFactors = F, header = F, sep = "")
Metadata <- read.delim("Metadata_Halibut_ENV_noNAnoWF.txt", header = T, stringsAsFactors = F)
ID <- FAM[1]
colnames(ID)[1] <- "ID"
Metadata_sorted <- inner_join(ID, Metadata)
#PCA CHR15
library(pcadapt)
Halib<- read.pcadapt("Halibut_CHR15.ped", type = "ped")
PCs <- pcadapt(Halib, K = 3, min.maf = 0.01)
PC_scores_pop <- as.data.frame(cbind(Metadata_sorted, PCs$scores))

plot(PCs, option = "manhattan" ) + theme_classic()
plot(PCs, option = "screeplot" ) + theme_classic()
#plot multiple axes
plot(PCs, option="scores", pop =  Metadata_sorted$Region, i = 1, j =2 ) + theme_classic() + geom_vline(xintercept = 0.00) + geom_vline(xintercept = 0.06)
plot(PCs, option="scores", pop =  Metadata_sorted$Region, i = 1, j =3 ) + theme_classic()
#Get inversion haplotypes
INV15_HAP1 <- Metadata_sorted$ID[which(PCs$scores[,1]<0.00)]
write.table(cbind(INV15_HAP1), "INV15_Hap1_inds.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(cbind(INV15_HAP1, INV15_HAP1), "INV15_Hap1_inds_plink.txt", sep = "\t", col.names = F, row.names = F, quote = F)

INV15_HAP2 <- Metadata_sorted$ID[which(PCs$scores[,1]>(0.06))]
write.table(cbind(INV15_HAP2), "INV15_Hap2_inds.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(cbind(INV15_HAP2, INV15_HAP2), "INV15_Hap2_inds_plink.txt", sep = "\t", col.names = F, row.names = F, quote = F)


INV15_HET <- Metadata_sorted$ID[which(PCs$scores[,1]<(0.06) & PCs$scores[,1]>0)]
write.table(cbind(INV15_HET), "INV15_HET_inds.txt", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(cbind(INV15_HET, INV15_HET), "INV15_HET_inds_plink.txt", sep = "\t", col.names = F, row.names = F, quote = F)
INV15_HET %in% INV15_HAP2

INV15_HAP1_Cluster <-  data.frame(cbind(INV15_HAP1, INV15_HAP1, rep("HAP1", length(INV15_HAP1))))
INV15_HAP2_Cluster <-  data.frame(cbind(INV15_HAP2, INV15_HAP2, rep("HAP2", length(INV15_HAP2))))
colnames(INV15_HAP1_Cluster) <- c("FID", "IID", "Cluster")
colnames(INV15_HAP2_Cluster) <- c("FID", "IID", "Cluster")
INV15_Cluster <- rbind(INV15_HAP1_Cluster, INV15_HAP2_Cluster)


#Estimate and plot FST
write.table(cbind(INV15_Cluster), "INV15_Cluster.txt", sep = "\t", col.names = F, row.names = F, quote = F)
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_CHR15 --fst --within INV15_Cluster.txt --out FST_CHR15_HAP1HAP2")
FST_CHR15_HAP1HAP2 <- read.delim("~/Desktop/Projects/Halibut_FinalRAD/FST_CHR15_HAP1HAP2.fst", stringsAsFactors=FALSE)
ggplot() + geom_point(data = FST_CHR15_HAP1HAP2, aes(POS, FST)) + theme_classic()  + geom_vline(xintercept = 5750000) + 
  geom_vline(xintercept = 12000000) + theme_classic()


#Heterozygosity by genotype
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV  --hardy --maf 0.01 --keep INV15_Hap1_inds_plink.txt --recode  --out CHR15_Hap1 --chr 15")
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV  --hardy --maf 0.01 --keep INV15_Hap2_inds_plink.txt --recode  --out CHR15_Hap2 --chr 15")
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV  --hardy --maf 0.01 --keep INV15_HET_inds_plink.txt --recode  --out CHR15_Het --chr 15")


INV15_Hap1_het  <- read.csv("~/Desktop/Projects/Halibut_FinalRAD/CHR15_Hap1.hwe", sep="")
INV15_Hap2_het  <- read.csv("~/Desktop/Projects/Halibut_FinalRAD/CHR15_Hap2.hwe", sep="")
INV15_Het_het  <- read.csv("~/Desktop/Projects/Halibut_FinalRAD/CHR15_Het.hwe", sep="")

INV15_map  <- read.csv("~/Desktop/Projects/Halibut_FinalRAD/CHR15.map", sep="", header = F)
colnames(INV15_map) <- c("Chrom", "SNP", "CM", "BP")
INV15_Hap1_het <- inner_join(INV15_Hap1_het, INV15_map)
INV15_Hap2_het <- inner_join(INV15_Hap2_het, INV15_map)
INV15_Het_het <- inner_join(INV15_Het_het, INV15_map)

ggplot() + geom_smooth(data = INV15_Hap1_het, aes(x = BP, y = O.HET., colour = "HAP1"), method = "loess",span = 0.03, se = F ) +
  geom_smooth(data = INV15_Hap2_het, aes(x = BP, y = O.HET., colour = "HAP2"), method = "loess",span = 0.03, se = F ) +
  geom_smooth(data = INV15_Het_het, aes(x = BP, y = O.HET., colour = "HET"), method = "loess",span = 0.03, se = F ) + geom_vline(xintercept = 5750000) + 
  geom_vline(xintercept = 12000000) + theme_classic()


#### compare heterozygosity by chromosome
CHR15_INV_Het_Het <- INV15_Het_het[which(INV15_Het_het$BP > 5750000 & INV15_Het_het$BP < 12000000),]
CHR15_NOTINV_NotINV_Het_Het <-INV15_Het_het[!(INV15_Het_het$SNP %in% CHR15_INV_Het_Het$SNP),]
wilcox.test(CHR15_INV_Het_Het$O.HET., CHR15_NOTINV_NotINV_Het_Het$O.HET.)
mean(CHR15_INV_Het_Het$O.HET.)
mean(CHR15_NOTINV_NotINV_Het_Het$O.HET.)


CHR15_INV_Hap2_Het <- INV15_Hap2_het[which(INV15_Hap2_het$BP > 5750000 & INV15_Hap2_het$BP < 12000000),]
CHR15_NOTINV_NotINV_Hap2_Het <-INV15_Hap2_het[!(INV15_Hap2_het$SNP %in% CHR15_INV_Hap2_Het$SNP),]
wilcox.test(CHR15_INV_Hap2_Het$O.HET., CHR15_NOTINV_NotINV_Hap2_Het$O.HET.)
mean(CHR15_INV_Hap2_Het$O.HET.)
mean(CHR15_NOTINV_NotINV_Hap2_Het$O.HET.)

CHR15_INV_Hap1_Het <- INV15_Hap1_het[which(INV15_Hap1_het$BP > 5750000 & INV15_Hap1_het$BP < 12000000),]
CHR15_NOTINV_NotINV_Het_Het <-INV15_Hap1_het[!(INV15_Hap1_het$SNP %in% CHR15_INV_Het_Het$SNP),]
wilcox.test(CHR15_INV_Het_Het$O.HET., CHR15_NOTINV_NotINV_Het_Het$O.HET.)
mean(CHR15_INV_Hap1_Het$O.HET.)
mean(CHR15_NOTINV_NotINV_Het_Het$O.HET.)

#LD chr15 comparisons ###
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV  --r2 square --maf 0.01 --keep INV15_Hap1_inds_plink.txt --recode  --out CHR15_Hap1_LD --chr 15")
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV  --r2 square --maf 0.01 --keep INV15_Hap2_inds_plink.txt --recode  --out CHR15_Hap2_LD --chr 15")
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV  --r2 square --maf 0.01 --keep INV15_HET_inds_plink.txt --recode  --out CHR15_Het_LD --chr 15")
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV  --r2 square --maf 0.01  --recode  --out CHR15_LD --chr 15")


#LD chr15 ###
library(data.table)

LD_matrix_CHR15 <- fread("CHR15_LD.ld", header = F)
map=data.table::fread("CHR15_LD.map")
colnames(LD_matrix_CHR15)=map$V2 #change column names to SNP names

#Set as matrix
LD_matrix_CHR15=as.matrix(LD_matrix_CHR15) 
#Change diagnoal to NAs
diag(LD_matrix_CHR15)<-NA


#Start and end of each matrix to calculate mean
#Run with full matrix  (make sure only diagonal is NA; ie. not missing upper or lower matrix values)

#Start at column 25 and get mean pairwise LD between the previous 25 and next 25 SNPs  (50 SNPs)

result<-NULL

#Run loop
for(i in 1:(nrow(LD_matrix_CHR15)-50)){
  row_end=i+50
  column_start=i+25
  result[i]<-mean(LD_matrix_CHR15[i:row_end,column_start], na.rm=T)
}
print(result)
str(LD_matrix_CHR15)
str(result) #Results will be shorter than LD matrix by 50 rows

y=nrow(LD_matrix_CHR15)-25 #Get number of rows in LD matrix wtihout last 25 (remove last 25 SNPs)
loci=map$V2[26:y] #Get SNP IDs starting at SNP 25 and going to end with last 25 removed (remove first 25 SNPs + last 25 SNPs)
all_LD_CHR15=cbind(result, loci) #Bind resutls from matrix calculations with Loci info

all_LD_CHR15 #Has Mean LD (for window) with SNP name

#Merge loci LD/name with Locus position for plot
data_for_plot=merge(all_LD_CHR15, map, by.x=2, by.y=2)
head(data_for_plot)

#Data for plot
data_for_plot=data_for_plot[,c(1,2,3,5)]
colnames(data_for_plot)=c("Locus", "Mean_LD", "Chr",  "Position")
#Change to numeric
data_for_plot$Position=as.numeric(as.character(data_for_plot$Position))
data_for_plot$Mean_LD=as.numeric(as.character(data_for_plot$Mean_LD))

str(data_for_plot)

#Quick plot
plot(data_for_plot$Position,data_for_plot$Mean_LD) #NOTE** see plot for nicer ggplot script BELOW
head(data_for_plot)

data_CHR15=data_for_plot

data_CHR15$SNP <- as.character(data_CHR15$Locus)


LD_FST <- inner_join(FST_CHR15_HAP1HAP2, data_CHR15, by = "SNP")
cor.test(LD_FST$FST, LD_FST$Mean_LD)
ggplot() + geom_point(data = LD_FST, aes(x = FST, y = Mean_LD))

LD_FST[which(LD_FST$POS > 5750000 & LD_FST$POS < 12000000),]
LD_FST_INV <- LD_FST[which(LD_FST$POS > 5750000 & LD_FST$POS < 12000000),]
mean(LD_FST_INV$FST, na.rm = T)
mean(LD_FST_INV$Mean_LD, na.rm = T)
cor.test(LD_FST_INV$FST, LD_FST_INV$Mean_LD)
LD_FST_OUTINV <- LD_FST[!(LD_FST$SNP %in% LD_FST_INV$SNP),]
mean(LD_FST_OUTINV$FST, na.rm = T)
mean(LD_FST_OUTINV$Mean_LD, na.rm = T)
cor.test(LD_FST_OUTINV$FST, LD_FST_OUTINV$Mean_LD)
wilcox.test(LD_FST_OUTINV$Mean_LD, LD_FST_INV$Mean_LD)
wilcox.test(LD_FST_OUTINV$FST, LD_FST_INV$FST)

#plot
ggplot() + geom_point(data = LD_FST_INV, aes(x = FST, y = Mean_LD, colour = "INV")) +
  geom_smooth(data = LD_FST_INV, aes(x = FST, y = Mean_LD, colour = "INV"),method='lm', formula= y~x ) +
  geom_point(data = LD_FST_OUTINV, aes(x = FST, y = Mean_LD, colour = "NOT_INV")) +
  theme_classic()


##heatmap of putative inversion region - 4MB to 13MB
system("~/Desktop/Software/plink_mac_20200219/plink --file Halibut_finalgenome_mindp15_24chrom_ming90_HWE_noWF_allENV  --r2 square --maf 0.01 --from-kb 4000 --to-kb 13000 --recode  --out CHR15_LD_Range --chr 15")
CHR15LDmat  <- as.matrix(fread("CHR15_LD_Range.ld", data.table = F))
diag(CHR15LDmat) <- NA
CHR15LDMat_Heatmap <- CHR15LDmat
CHR15LDMat_Heatmap[upper.tri(CHR15LDMat_Heatmap, diag = T)] <- 0
heatpal <-c("gray95", "blue", "purple", "red2", "red3")
superheat(CHR15LDMat_Heatmap, heat.pal = heatpal, heat.pal.values = c(0, 0.25, 0.5, 0.75, 1))


#repeat 50-SNP LD calculation loop for HET, HAP1, HAP2 (just change input and output names) then plot

ggplot() + geom_smooth(data = data_CHR15, aes(x = Position, y = Mean_LD, colour = "All"), method = "loess", span = 0.05, se = F) +
  geom_smooth(data = data_CHR15_HAP1, aes(x = Position, y = Mean_LD, colour = "Hap1"), method = "loess", span = 0.05, se = F) +
  geom_smooth(data = data_CHR15_HAP2, aes(x = Position, y = Mean_LD, colour = "Hap2"), method = "loess", span = 0.05, se = F) +
  geom_smooth(data = data_CHR15_HET, aes(x = Position, y = Mean_LD, colour = "Het"), method = "loess", span = 0.05, se = F)+ geom_vline(xintercept = 5750000) + 
  geom_vline(xintercept = 12000000) + theme_classic()

CHR15_INV_HAP1_LD <- data_CHR15_HAP1[which(data_CHR15_HAP1$Position > 5750000 & data_CHR15_HAP1$Position < 12000000),]
CHR15_NOTINV_HAP1_LD <- data_CHR15_HAP1[!(data_CHR15_HAP1$SNP %in% CHR15_INV_HAP1_LD$SNP),]
mean(CHR15_INV_HAP1_LD$Mean_LD)
mean(CHR15_NOTINV_HAP1_LD$Mean_LD)
wilcox.test(CHR15_INV_HAP1_LD$Mean_LD, CHR15_NOTINV_HAP1_LD$Mean_LD)

CHR15_INV_HAP2_LD <- data_CHR15_HAP2[which(data_CHR15_HAP2$Position > 5750000 & data_CHR15_HAP2$Position < 12000000),]
CHR15_NOTINV_HAP2_LD <- data_CHR15_HAP2[!(data_CHR15_HAP2$SNP %in% CHR15_INV_HAP2_LD$SNP),]
mean(CHR15_INV_HAP2_LD$Mean_LD)
mean(CHR15_NOTINV_HAP2_LD$Mean_LD)
wilcox.test(CHR15_INV_HAP2_LD$Mean_LD, CHR15_NOTINV_HAP2_LD$Mean_LD)

CHR15_INV_HET_LD <- data_CHR15_HET[which(data_CHR15_HET$Position > 5750000 & data_CHR15_HET$Position < 12000000),]
CHR15_NOTINV_HET_LD <- data_CHR15_HET[!(data_CHR15_HET$SNP %in% CHR15_INV_HET_LD$SNP),]
mean(CHR15_INV_HET_LD$Mean_LD)
mean(CHR15_NOTINV_HET_LD$Mean_LD)
wilcox.test(CHR15_INV_HET_LD$Mean_LD, CHR15_NOTINV_HET_LD$Mean_LD)

CHR15_INV_LD <- data_CHR15[which(data_CHR15$Position > 5750000 & data_CHR15$Position < 12000000),]
CHR15_NOTINV_LD <- data_CHR15[!(data_CHR15$SNP %in% CHR15_INV_LD$SNP),]
mean(CHR15_INV_LD$Mean_LD)
mean(CHR15_NOTINV_LD$Mean_LD)
wilcox.test(CHR15_INV_LD$Mean_LD, CHR15_NOTINV_LD$Mean_LD)

#GO prep
system("bedtools intersect -b INV15_minimapaligned.bed -a GCF_009819705.1_fHipHip1.pri_genomic.gff")
INV15_region <- read.delim("~/Desktop/Projects/Halibut_FinalRAD/INV15_region.gff", header=FALSE)
INV15_region_genes_only <- INV15_region[INV15_region$V3 %in% "gene",]
#Keep only genes with symbols 
INV15_region_genes_only <- INV15_region_genes_only[!(grepl(pattern = "ID=gene-LOC1.*", INV15_region_genes_only$V9)),]
write.table(INV15_region_genes_only, "INV15_gene_only.gff", col.names = F, sep = "\t", row.names = F, quote = F)

#Map of putative inversion distribution
#it pie time
library(maps)
library(mapdata) # all your basemaps are here
library(mapplots) # for add.pie
library(gplots) # for colour range
library(rworldmap)
library(maptools)
library(lattice)
INV15_HAP1 <- Metadata_sorted[which(PCs$scores[,1]<0.00),]
INV15_HAP1$Hap <- "HAP1"

INV15_HAP2 <- Metadata_sorted[which(PCs$scores[,1]>(0.06)),]
INV15_HAP2$Hap <- "HAP2"

INV15_HET <- Metadata_sorted[which(PCs$scores[,1]<(0.06) & PCs$scores[,1]>0),]
INV15_HET$Hap <- "HET" 

Metadata_haplotypes <- rbind(INV15_HET, INV15_HAP2, INV15_HAP1)

Mean_lat <- Metadata_haplotypes %>% 
  group_by(Region ) %>% 
  summarise(mean_lat = mean(Lat))

Mean_lon <- Metadata_haplotypes %>% 
  group_by(Region ) %>% 
  summarise(mean_lon = mean(Lon))

Mean_coords <- inner_join(Mean_lat, Mean_lon)

Group_INVfreq <- data.frame(cbind(Metadata_haplotypes$Region, Metadata_haplotypes$Hap), stringsAsFactors = F)
INV_Freq_Table <- as.data.frame.matrix(prop.table(table(Group_INVfreq), 1) * 100)
INV_Freq_Table$Region <- rownames(INV_Freq_Table)
INV_Freq_Coords<- data.frame(inner_join(Mean_coords, INV_Freq_Table), stringsAsFactors = F)

#get map
Sample.Lat.lim=c(40,51)
Sample.Long.lim=c(-71,-43)
map("worldHires", xlim=Sample.Long.lim, ylim=Sample.Lat.lim, col="grey80",  fill=T, lwd=0.01,  resolution=0);map.axes()

#add pies
for(i in 1:nrow(INV_Freq_Coords)){
  add.pie(as.integer(INV_Freq_Coords[i,c("HET", "HAP1", "HAP2")]),
          x=INV_Freq_Coords$mean_lon[i],y=INV_Freq_Coords$mean_lat[i],labels="",border = F, radius = 0.6,lty = NULL, density = NULL, col=c("blue","orange", "red"))}

