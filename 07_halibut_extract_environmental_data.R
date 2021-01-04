Metadata <- read.delim("~/Desktop/Halib_allgenomealigned/Halibut_metadata_noNAcoords.txt", header = T, stringsAsFactors = F)

library(sdmpredictors)
library(wesanderson)
datasets <- list_datasets(terrestrial = FALSE, marine = TRUE)
datasets
layers <- list_layers(datasets)

ne.atlantic.ext <- extent(-72.5, 15, 41, 70) 

temp.max.bottom <- load_layers("BO2_tempmax_bdmax")
temp.max.bottom.crop <- crop(temp.max.bottom, ne.atlantic.ext)

temp.mean.bottom <- load_layers("BO2_tempmean_bdmax")
temp.mean.bottom.crop <- crop(temp.mean.bottom, ne.atlantic.ext)

temp.min.bottom <- load_layers("BO2_tempmin_bdmax")
temp.min.bottom.crop <- crop(temp.min.bottom, ne.atlantic.ext)


light.max.bottom <- load_layers("BO2_lightbotmax_bdmax")
light.max.bottom.crop <- crop(light.max.bottom, ne.atlantic.ext)

light.mean.bottom <- load_layers("BO2_lightbotmean_bdmax")
light.mean.bottom.crop <- crop(light.mean.bottom, ne.atlantic.ext)

light.min.bottom <- load_layers("BO2_lightbotmin_bdmax")
light.min.bottom.crop <- crop(light.min.bottom, ne.atlantic.ext)

o2.max.bottom <- load_layers("BO2_dissoxmax_bdmax")
o2.max.bottom.crop <- crop(o2.max.bottom, ne.atlantic.ext)

o2.mean.bottom <- load_layers("BO2_dissoxmean_bdmax")
o2.mean.bottom.crop <- crop(o2.mean.bottom, ne.atlantic.ext)

o2.min.bottom <- load_layers("BO2_dissoxmin_bdmax")
o2.min.bottom.crop <- crop(o2.min.bottom, ne.atlantic.ext)

salt.max.bottom <- load_layers("BO2_salinitymax_bdmax")
salt.max.bottom.crop <- crop(salt.max.bottom, ne.atlantic.ext)

salt.mean.bottom <- load_layers("BO2_salinitymean_bdmax")
salt.mean.bottom.crop <- crop(salt.mean.bottom, ne.atlantic.ext)

salt.min.bottom <- load_layers("BO2_salinitymin_bdmax")
salt.min.bottom.crop <- crop(salt.min.bottom, ne.atlantic.ext)


my.colors = colorRampPalette(wes_palette("Zissou1")) 
plot(temp.max.bottom.crop,col=my.colors(1000), box=FALSE) 
plot(temp.mean.bottom.crop,col=my.colors(1000), box=FALSE) 
plot(temp.min.bottom.crop,col=my.colors(1000), box=FALSE) 
plot(light.max.bottom.crop, col=my.colors(1000), box=FALSE) 
plot(light.mean.bottom.crop, col=my.colors(1000), box=FALSE) 
plot(light.min.bottom.crop, col=my.colors(1000), box=FALSE) 
plot(o2.max.bottom.crop,col=my.colors(1000), box=FALSE) 
plot(o2.mean.bottom.crop,col=my.colors(1000), box=FALSE) 
plot(o2.min.bottom.crop,col=my.colors(1000), box=FALSE) 
plot(salt.max.bottom.crop,col=my.colors(1000), box=FALSE) 
plot(salt.mean.bottom.crop,col=my.colors(1000), box=FALSE) 
plot(salt.min.bottom.crop,col=my.colors(1000), box=FALSE) 

environment.bottom <- load_layers(c("BO2_tempmax_bdmean", "BO2_tempmean_bdmean", "BO2_tempmin_bdmean", "BO2_lightbotmax_bdmean", "BO2_lightbotmean_bdmean", "BO2_lightbotmin_bdmean", "BO2_dissoxmax_bdmean", "BO2_dissoxmean_bdmean", "BO2_dissoxmin_bdmean", "BO2_salinitymax_bdmean", "BO2_salinitymean_bdmean", "BO2_salinitymin_bdmean"))


bathymetry <- load_layers("BO_bathymean") 
bathymetry[,2:3]
my.sites_halb <- data.frame(cbind(ID = Metadata$ID ,Lon= as.numeric(Metadata$Lon), Lat = as.numeric(Metadata$Lat))) 
my.sites_halb$Lat <- as.numeric(as.character(my.sites_halb$Lat))
my.sites_halb$Lon <- as.numeric(as.character(my.sites_halb$Lon))


my.sites.environment <- data.frame(ID=my.sites_halb$ID , depth=extract(bathymetry,my.sites_halb[,2:3]) , extract(environment.bottom,my.sites_halb[,2:3]) ) 
my.sites.environment 

library(dplyr)
Halibut_Metadata_ENV <- inner_join(Metadata, my.sites.environment)
#Remove those exhibiting high structuring prior to standardize! 
Halibut_Metadata_ENV <- Halibut_Metadata_ENV[!(Halibut_Metadata_ENV$ID %in% High_structure_halibuts$V1),]

Halibut_Metadata_ENV$BO2_tempmax_bdmean_sd <- Halibut_Metadata_ENV$BO2_tempmax_bdmean/sd(Halibut_Metadata_ENV$BO2_tempmax_bdmean)
Halibut_Metadata_ENV$BO2_tempmean_bdmean_sd <- Halibut_Metadata_ENV$BO2_tempmean_bdmean/sd(Halibut_Metadata_ENV$BO2_tempmean_bdmean)
Halibut_Metadata_ENV$BO2_tempmin_bdmean_sd <- Halibut_Metadata_ENV$BO2_tempmin_bdmean/sd(Halibut_Metadata_ENV$BO2_tempmin_bdmean)
Halibut_Metadata_ENV$BO2_dissoxmax_bdmean_sd <- Halibut_Metadata_ENV$BO2_dissoxmax_bdmean/sd(Halibut_Metadata_ENV$BO2_dissoxmax_bdmean)
Halibut_Metadata_ENV$BO2_dissoxmean_bdmean_sd <- Halibut_Metadata_ENV$BO2_dissoxmean_bdmean/sd(Halibut_Metadata_ENV$BO2_dissoxmean_bdmean)
Halibut_Metadata_ENV$BO2_dissoxmin_bdmean_sd <- Halibut_Metadata_ENV$BO2_dissoxmin_bdmean/sd(Halibut_Metadata_ENV$BO2_dissoxmin_bdmean)
Halibut_Metadata_ENV$BO2_salinitymax_bdmean_sd <- Halibut_Metadata_ENV$BO2_salinitymax_bdmean/sd(Halibut_Metadata_ENV$BO2_salinitymax_bdmean)
Halibut_Metadata_ENV$BO2_salinitymean_bdmean_sd <- Halibut_Metadata_ENV$BO2_salinitymean_bdmean/sd(Halibut_Metadata_ENV$BO2_salinitymean_bdmean)
Halibut_Metadata_ENV$BO2_salinitymin_bdmean_sd <- Halibut_Metadata_ENV$BO2_salinitymin_bdmean/sd(Halibut_Metadata_ENV$BO2_salinitymin_bdmean)


##PCA
TEMP_PC <- prcomp(cbind(Halibut_Metadata_ENV$BO2_tempmax_bdmean_sd, Halibut_Metadata_ENV$BO2_tempmean_bdmean_sd, Halibut_Metadata_ENV$BO2_tempmin_bdmean_sd))
class(TEMP_mat)
TEMP_PC <- prcomp(TEMP_mat)
TEMP_PC1 <- TEMP_PC$x[,1]

O2_PC <- prcomp(cbind(Halibut_Metadata_ENV$BO2_dissoxmax_bdmean_sd, Halibut_Metadata_ENV$BO2_dissoxmean_bdmean_sd, Halibut_Metadata_ENV$BO2_dissoxmin_sd))
O2_PC1 <- O2_PC$x[,1]

SALT_PC <- prcomp(cbind(Halibut_Metadata_ENV$BO2_salinitymax_bdmean_sd, Halibut_Metadata_ENV$BO2_salinitymean_bdmean_sd, Halibut_Metadata_ENV$BO2_salinitymin_sd))
SALT_PC1 <- SALT_PC$x[,1]
Halibut_Metadata_ENV$TEMP_PC <- TEMP_PC1
Halibut_Metadata_ENV$O2_PC <- O2_PC1
Halibut_Metadata_ENV$SALT_PC <- SALT_PC1
write.table(x = Halibut_Metadata_ENV, file = "Halibut_Metadata_ENV_Complete.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
