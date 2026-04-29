setwd("/blue/ferrallm/antonihalka")
###Preprocessing of Pancreas samples
###Samples from GSE148073, GSE154778, GSE154778

#Load libraries
library('ggplot2')
library('Seurat')
library('patchwork')
library('forcats')
library('tidyverse')
library('data.table')
library('glmGamPoi')
library('future')
library('harmony')
library('clustree')
library('plyr')
library('hdf5r')
library('SingleR')
library(dplyr)
library(BiocParallel)
#library(ArrayExpress) #Just in case we need it later
library('celldex')
library(clusterProfiler)
library(org.Hs.eg.db)
#library(enrichplot)
library(msigdbr)
library(ggrepel)

# In R script or R console
install.packages("Seurat")
install.packages("SeuratObject")

library(Seurat)
library(SeuratObject)

packageVersion("Seurat")
packageVersion("SeuratObject")

# Increase max size for globals to 4 GB then 20GB
BiocParallel::register(BiocParallel::MulticoreParam(workers = 8)) 
options(future.globals.maxSize = 50 * 1024^3)

## Set up Directories
dir <- "/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/rds/"
figdir <- "/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/"

###GSE205013#####################################

###P01_Bam
P01file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204109/P01_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P01.counts <- Read10X(data.dir=P01file)

# create a Seurat object based on the scRNA-seq data
P01 <- CreateSeuratObject(counts=P01.counts, project="P01_GSM6204109")
#add metadata
P01@meta.data[["stage"]] <- "4"
P01@meta.data[["tissue type"]] <- "Metastatsis"
P01@meta.data[["site"]] <- "Liver"
P01@meta.data[["treatment"]] <- "Untreated"
P01@meta.data[["treatment.type"]] <- "Untreated"
P01@meta.data[["sex"]] <- "M"
P01@meta.data[["age"]] <- "78"

#save rds file
saveRDS(P01, paste(dir,"/P01_GSM6204109.rds",sep=""))

###P02_Bam
P02file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204110/P02_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P02.counts <- Read10X(data.dir=P02file)

# create a Seurat object based on the scRNA-seq data
P02 <- CreateSeuratObject(counts=P02.counts, project="P02_GSM6204110")
#add metadata
P02@meta.data[["stage"]] <- "4"
P02@meta.data[["tissue type"]] <- "Metastatsis"
P02@meta.data[["site"]] <- "Liver"
P02@meta.data[["treatment"]] <- "Untreated"
P02@meta.data[["treatment.type"]] <- "Untreated"
P02@meta.data[["sex"]] <- "M"
P02@meta.data[["age"]] <- "74"

#save rds file
saveRDS(P02, paste(dir,"/P02_GSM6204110.rds",sep=""))


###P03_Bam
P03file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204111/P03_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P03.counts <- Read10X(data.dir=P03file)

P03 <- CreateSeuratObject(counts=P03.counts, project="P03_GSM6204111")
P03@meta.data[["stage"]] <- "4"
P03@meta.data[["tissue type"]] <- "Primary"
P03@meta.data[["site"]] <- "Pancreas"
P03@meta.data[["treatment"]] <- "Treated"
P03@meta.data[["treatment.type"]] <- "chemo"
P03@meta.data[["sex"]] <- "M"
P03@meta.data[["age"]] <- "81"

saveRDS(P03, paste(dir,"/P03_GSM6204111.rds",sep=""))

###P04_Bam
P04file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204112/P04_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P04.counts <- Read10X(data.dir=P04file)

P04 <- CreateSeuratObject(counts=P04.counts, project="P04_GSM6204112")
P04@meta.data[["stage"]] <- "1"
P04@meta.data[["tissue type"]] <- "Primary"
P04@meta.data[["site"]] <- "Pancreas"
P04@meta.data[["treatment"]] <- "Untreated"
P04@meta.data[["treatment.type"]] <- "Untreated"
P04@meta.data[["sex"]] <- "F"
P04@meta.data[["age"]] <- "85"

saveRDS(P04, paste(dir,"/P04_GSM6204112.rds",sep=""))

###P05_Bam
P05file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204113/P05_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P05.counts <- Read10X(data.dir=P05file)

P05 <- CreateSeuratObject(counts=P05.counts, project="P05_GSM6204113")
P05@meta.data[["stage"]] <- "3"
P05@meta.data[["tissue type"]] <- "Primary"
P05@meta.data[["site"]] <- "Pancreas"
P05@meta.data[["treatment"]] <- "Untreated"
P05@meta.data[["treatment.type"]] <- "Untreated"
P05@meta.data[["sex"]] <- "F"
P05@meta.data[["age"]] <- "69"

saveRDS(P05, paste(dir,"/P05_GSM6204113.rds",sep=""))

###P06_Bam
P06file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204114/P06_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P06.counts <- Read10X(data.dir=P06file)

P06 <- CreateSeuratObject(counts=P06.counts, project="P06_GSM6204114")
P06@meta.data[["stage"]] <- "3"
P06@meta.data[["tissue type"]] <- "Primary"
P06@meta.data[["site"]] <- "Pancreas"
P06@meta.data[["treatment"]] <- "Treated"
P06@meta.data[["treatment.type"]] <- "chemo"
P06@meta.data[["sex"]] <- "M"
P06@meta.data[["age"]] <- "55"

saveRDS(P06, paste(dir,"/P06_GSM6204114.rds",sep=""))

###P07_Bam
P07file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204115/P07_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P07.counts <- Read10X(data.dir=P07file)

P07 <- CreateSeuratObject(counts=P07.counts, project="P07_GSM6204115")
P07@meta.data[["stage"]] <- "1"
P07@meta.data[["tissue type"]] <- "Primary"
P07@meta.data[["site"]] <- "Pancreas"
P07@meta.data[["treatment"]] <- "Untreated"
P07@meta.data[["treatment.type"]] <- "Untreated"
P07@meta.data[["sex"]] <- "F"
P07@meta.data[["age"]] <- "70"

saveRDS(P07, paste(dir,"/P07_GSM6204115.rds",sep=""))

###P08_Bam
P08file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204116/P08_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P08.counts <- Read10X(data.dir=P08file)

P08 <- CreateSeuratObject(counts=P08.counts, project="P08_GSM6204116")
P08@meta.data[["stage"]] <- "2"
P08@meta.data[["tissue type"]] <- "Primary"
P08@meta.data[["site"]] <- "Pancreas"
P08@meta.data[["treatment"]] <- "Treated"
P08@meta.data[["treatment.type"]] <- "chemo"
P08@meta.data[["sex"]] <- "F"
P08@meta.data[["age"]] <- "87"

saveRDS(P08, paste(dir,"/P08_GSM6204116.rds",sep=""))

###P09_Bam
P09file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204117/P09_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P09.counts <- Read10X(data.dir=P09file)

P09 <- CreateSeuratObject(counts=P09.counts, project="P09_GSM6204117")
P09@meta.data[["stage"]] <- "4"
P09@meta.data[["tissue type"]] <- "Primary"
P09@meta.data[["site"]] <- "Pancreas"
P09@meta.data[["treatment"]] <- "Untreated"
P09@meta.data[["treatment.type"]] <- "Untreated"
P09@meta.data[["sex"]] <- "M"
P09@meta.data[["age"]] <- "63"

saveRDS(P09, paste(dir,"/P09_GSM6204117.rds",sep=""))

###P10_Bam
P10file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204118/P10_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P10.counts <- Read10X(data.dir=P10file)

P10 <- CreateSeuratObject(counts=P10.counts, project="P10_GSM6204118")
P10@meta.data[["stage"]] <- "2"
P10@meta.data[["tissue type"]] <- "Primary"
P10@meta.data[["site"]] <- "Pancreas"
P10@meta.data[["treatment"]] <- "Treated"
P10@meta.data[["treatment.type"]] <- "chemo"
P10@meta.data[["sex"]] <- "F"
P10@meta.data[["age"]] <- "48"

saveRDS(P10, paste(dir,"/P10_GSM6204118.rds",sep=""))

###P11_Bam
P11file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204119/P11_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P11.counts <- Read10X(data.dir=P11file)

# create a Seurat object based on the scRNA-seq data
P11 <- CreateSeuratObject(counts=P11.counts, project="P11_GSM6204119")
#add metadata
P11@meta.data[["stage"]] <- "4"
P11@meta.data[["tissue type"]] <- "Metastasis"
P11@meta.data[["site"]] <- "Liver"
P11@meta.data[["treatment"]] <- "Untreated"
P11@meta.data[["treatment.type"]] <- "Untreated"
P11@meta.data[["sex"]] <- "F"
P11@meta.data[["age"]] <- "73"

#save rds file
saveRDS(P11, paste(dir,"/P11_GSM6204119.rds",sep=""))

###P12_Bam
P12file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204120/P12_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P12.counts <- Read10X(data.dir=P12file)

# create a Seurat object based on the scRNA-seq data
P12 <- CreateSeuratObject(counts=P12.counts, project="P12_GSM6204120")
#add metadata
P12@meta.data[["stage"]] <- "4"
P12@meta.data[["tissue type"]] <- "Primary"
P12@meta.data[["site"]] <- "Pancreas"
P12@meta.data[["treatment"]] <- "Treated"
P12@meta.data[["treatment.type"]] <- "chemo"
P12@meta.data[["sex"]] <- "M"
P12@meta.data[["age"]] <- "68"

#save rds file
saveRDS(P12, paste(dir,"/P12_GSM6204120.rds",sep=""))

###P13_Bam
P13file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204121/P13_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P13.counts <- Read10X(data.dir=P13file)

# create a Seurat object based on the scRNA-seq data
P13 <- CreateSeuratObject(counts=P13.counts, project="P13_GSM6204121")
#add metadata
P13@meta.data[["stage"]] <- "3"
P13@meta.data[["tissue type"]] <- "Primary"
P13@meta.data[["site"]] <- "Pancreas"
P13@meta.data[["treatment"]] <- "Untreated"
P13@meta.data[["treatment.type"]] <- "Untreated"
P13@meta.data[["sex"]] <- "M"
P13@meta.data[["age"]] <- "78"

#save rds file
saveRDS(P13, paste(dir,"/P13_GSM6204121.rds",sep=""))

###P14_Bam
P14file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204122/P14_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P14.counts <- Read10X(data.dir=P14file)

# create a Seurat object based on the scRNA-seq data
P14 <- CreateSeuratObject(counts=P14.counts, project="P14_GSM6204122")
#add metadata
P14@meta.data[["stage"]] <- "1"
P14@meta.data[["tissue type"]] <- "Primary"
P14@meta.data[["site"]] <- "Pancreas"
P14@meta.data[["treatment"]] <- "Treated"
P14@meta.data[["treatment.type"]] <- "chemo"
P14@meta.data[["sex"]] <- "F"
P14@meta.data[["age"]] <- "67"

#save rds file
saveRDS(P14, paste(dir,"/P14_GSM6204122.rds",sep=""))

###P15_Bam
P15file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204123/P15_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P15.counts <- Read10X(data.dir=P15file)

# create a Seurat object based on the scRNA-seq data
P15 <- CreateSeuratObject(counts=P15.counts, project="P15_GSM6204123")
#add metadata
P15@meta.data[["stage"]] <- "3"
P15@meta.data[["tissue type"]] <- "Primary"
P15@meta.data[["site"]] <- "Pancreas"
P15@meta.data[["treatment"]] <- "Untreated"
P15@meta.data[["treatment.type"]] <- "Untreated"
P15@meta.data[["sex"]] <- "F"
P15@meta.data[["age"]] <- "58"

#save rds file
saveRDS(P15, paste(dir,"/P15_GSM6204123.rds",sep=""))

###P16_Bam
P16file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204124/P16_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P16.counts <- Read10X(data.dir=P16file)

# create a Seurat object based on the scRNA-seq data
P16 <- CreateSeuratObject(counts=P16.counts, project="P16_GSM6204124")
#add metadata
P16@meta.data[["stage"]] <- "4"
P16@meta.data[["tissue type"]] <- "Metastasis"
P16@meta.data[["site"]] <- "Liver"
P16@meta.data[["treatment"]] <- "Untreated"
P16@meta.data[["treatment.type"]] <- "Untreated"
P16@meta.data[["sex"]] <- "F"
P16@meta.data[["age"]] <- "69"

#save rds file
saveRDS(P16, paste(dir,"/P16_GSM6204124.rds",sep=""))

###P17_Bam
P17file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204125/P17_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P17.counts <- Read10X(data.dir=P17file)

# create a Seurat object based on the scRNA-seq data
P17 <- CreateSeuratObject(counts=P17.counts, project="P17_GSM6204125")
#add metadata
P17@meta.data[["stage"]] <- "4"
P17@meta.data[["tissue type"]] <- "Metastasis"
P17@meta.data[["site"]] <- "Liver"
P17@meta.data[["treatment"]] <- "Treated"
P17@meta.data[["treatment.type"]] <- "chemo"
P17@meta.data[["sex"]] <- "M"
P17@meta.data[["age"]] <- "64"

#save rds file
saveRDS(P17, paste(dir,"/P17_GSM6204125.rds",sep=""))

###P19_Bam
P19file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204127/P19_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P19.counts <- Read10X(data.dir=P19file)

# create a Seurat object based on the scRNA-seq data
P19 <- CreateSeuratObject(counts=P19.counts, project="P19_GSM6204127")
#add metadata
P19@meta.data[["stage"]] <- "1"
P19@meta.data[["tissue type"]] <- "Primary"
P19@meta.data[["site"]] <- "Pancreas"
P19@meta.data[["treatment"]] <- "Untreated"
P19@meta.data[["treatment.type"]] <- "Untreated"
P19@meta.data[["sex"]] <- "F"
P19@meta.data[["age"]] <- "60"

#save rds file
saveRDS(P19, paste(dir,"/P19_GSM6204127.rds",sep=""))

###P20_Bam
P20file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204128/P20_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P20.counts <- Read10X(data.dir=P20file)

# create a Seurat object based on the scRNA-seq data
P20 <- CreateSeuratObject(counts=P20.counts, project="P20_GSM6204128")
#add metadata
P20@meta.data[["stage"]] <- "4"
P20@meta.data[["tissue type"]] <- "Primary"
P20@meta.data[["site"]] <- "Pancreas"
P20@meta.data[["treatment"]] <- "Untreated"
P20@meta.data[["treatment.type"]] <- "Untreated"
P20@meta.data[["sex"]] <- "F"
P20@meta.data[["age"]] <- "66"

#save rds file
saveRDS(P20, paste(dir,"/P20_GSM6204128.rds",sep=""))

###P21_Bam
P21file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204129/P21_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P21.counts <- Read10X(data.dir=P21file)

# create a Seurat object based on the scRNA-seq data
P21 <- CreateSeuratObject(counts=P21.counts, project="P21_GSM6204129")
#add metadata
P21@meta.data[["stage"]] <- "4"
P21@meta.data[["tissue type"]] <- "Metastasis"
P21@meta.data[["site"]] <- "Liver"
P21@meta.data[["treatment"]] <- "Untreated"
P21@meta.data[["treatment.type"]] <- "Untreated"
P21@meta.data[["sex"]] <- "M"
P21@meta.data[["age"]] <- "69"

#save rds file
saveRDS(P21, paste(dir,"/P21_GSM6204129.rds",sep=""))

###P24_Bam
P24file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204132/P24_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P24.counts <- Read10X(data.dir=P24file)

# create a Seurat object based on the scRNA-seq data
P24 <- CreateSeuratObject(counts=P24.counts, project="P24_GSM6204132")
#add metadata
P24@meta.data[["stage"]] <- "4"
P24@meta.data[["tissue type"]] <- "Metastasis"
P24@meta.data[["site"]] <- "Liver"
P24@meta.data[["treatment"]] <- "Untreated"
P24@meta.data[["treatment.type"]] <- "Untreated"
P24@meta.data[["sex"]] <- "F"
P24@meta.data[["age"]] <- "65"

#save rds file
saveRDS(P24, paste(dir,"/P24_GSM6204132.rds",sep=""))

###P25_Bam
P25file <- "/blue/ferrallm/antonihalka/Pancreas/Werba_et_al/GSM6204133/P25_bam/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P25.counts <- Read10X(data.dir=P25file)

# create a Seurat object based on the scRNA-seq data
P25 <- CreateSeuratObject(counts=P25.counts, project="P25_GSM6204133")
#add metadata
P25@meta.data[["stage"]] <- "4"
P25@meta.data[["tissue type"]] <- "Metastasis"
P25@meta.data[["site"]] <- "Liver"
P25@meta.data[["treatment"]] <- "Untreated"
P25@meta.data[["treatment.type"]] <- "Untreated"
P25@meta.data[["sex"]] <- "M"
P25@meta.data[["age"]] <- "48"

#save rds file
saveRDS(P25, paste(dir,"/P25_GSM6204133.rds",sep=""))

###GSE154778#####################################

###P03_Lin
P03_Linfile <- "/blue/ferrallm/antonihalka/Pancreas/Lin/P03/SRX8778047/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P03_Lin.counts <- Read10X(data.dir=P03_Linfile)

# create a Seurat object based on the scRNA-seq data
P03_Lin <- CreateSeuratObject(counts=P03_Lin.counts, project="P03_Lin")
#add metadata
P03_Lin@meta.data[["stage"]] <- "2"
P03_Lin@meta.data[["tissue type"]] <- "Primary"
P03_Lin@meta.data[["site"]] <- "Pancreas"
P03_Lin@meta.data[["treatment"]] <- "Unspecified"
P03_Lin@meta.data[["treatment.type"]] <- "Unspecified"
P03_Lin@meta.data[["sex"]] <- "M"
P03_Lin@meta.data[["age"]] <- "72"

#save rds file
saveRDS(P03_Lin, paste(dir,"/P03_Lin.rds",sep=""))

###P04_Lin
P04_Linfile <- "/blue/ferrallm/antonihalka/Pancreas/Lin/P04/SRX8778048/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P04_Lin.counts <- Read10X(data.dir=P04_Linfile)

# create a Seurat object based on the scRNA-seq data
P04_Lin <- CreateSeuratObject(counts=P04_Lin.counts, project="P04_Lin")
#add metadata
P04_Lin@meta.data[["stage"]] <- "2"
P04_Lin@meta.data[["tissue type"]] <- "Primary"
P04_Lin@meta.data[["site"]] <- "Pancreas"
P04_Lin@meta.data[["treatment"]] <- "Unspecified"
P04_Lin@meta.data[["treatment.type"]] <- "Unspecified"
P04_Lin@meta.data[["sex"]] <- "M"
P04_Lin@meta.data[["age"]] <- "57"

#save rds file
saveRDS(P04_Lin, paste(dir,"/P04_Lin.rds",sep=""))

###P05_Lin
P05_Linfile <- "/blue/ferrallm/antonihalka/Pancreas/Lin/P05/SRX8778049/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P05_Lin.counts <- Read10X(data.dir=P05_Linfile)

P05_Lin <- CreateSeuratObject(counts=P05_Lin.counts, project="P05_Lin")
P05_Lin@meta.data[["stage"]] <- "2"
P05_Lin@meta.data[["tissue type"]] <- "Primary"
P05_Lin@meta.data[["site"]] <- "Pancreas"
P05_Lin@meta.data[["treatment"]] <- "Unspecified"
P05_Lin@meta.data[["treatment.type"]] <- "Unspecified"
P05_Lin@meta.data[["sex"]] <- "M"
P05_Lin@meta.data[["age"]] <- "80"

saveRDS(P05_Lin, paste(dir,"/P05_Lin.rds",sep=""))


###P06_Lin
P06_Linfile <- "/blue/ferrallm/antonihalka/Pancreas/Lin/P06/SRX8778050/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P06_Lin.counts <- Read10X(data.dir=P06_Linfile)

P06_Lin <- CreateSeuratObject(counts=P06_Lin.counts, project="P06_Lin")
P06_Lin@meta.data[["stage"]] <- "2"
P06_Lin@meta.data[["tissue type"]] <- "Primary"
P06_Lin@meta.data[["site"]] <- "Pancreas"
P06_Lin@meta.data[["treatment"]] <- "Unspecified"
P06_Lin@meta.data[["treatment.type"]] <- "Unspecified"
P06_Lin@meta.data[["sex"]] <- "M"
P06_Lin@meta.data[["age"]] <- "41"

saveRDS(P06_Lin, paste(dir,"/P06_Lin.rds",sep=""))


###P07_Lin
P07_Linfile <- "/blue/ferrallm/antonihalka/Pancreas/Lin/P07/SRX8778051/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P07_Lin.counts <- Read10X(data.dir=P07_Linfile)

P07_Lin <- CreateSeuratObject(counts=P07_Lin.counts, project="P07_Lin")
P07_Lin@meta.data[["stage"]] <- "2"
P07_Lin@meta.data[["tissue type"]] <- "Primary"
P07_Lin@meta.data[["site"]] <- "Pancreas"
P07_Lin@meta.data[["treatment"]] <- "Unspecified"
P07_Lin@meta.data[["treatment.type"]] <- "Unspecified"
P07_Lin@meta.data[["sex"]] <- "F"
P07_Lin@meta.data[["age"]] <- "56"

saveRDS(P07_Lin, paste(dir,"/P07_Lin.rds",sep=""))


###P08_Lin
P08_Linfile <- "/blue/ferrallm/antonihalka/Pancreas/Lin/P08/SRX8778052/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
P08_Lin.counts <- Read10X(data.dir=P08_Linfile)

P08_Lin <- CreateSeuratObject(counts=P08_Lin.counts, project="P08_Lin")
P08_Lin@meta.data[["stage"]] <- "2"
P08_Lin@meta.data[["tissue type"]] <- "Primary"
P08_Lin@meta.data[["site"]] <- "Pancreas"
P08_Lin@meta.data[["treatment"]] <- "Unspecified"
P08_Lin@meta.data[["treatment.type"]] <- "Unspecified"
P08_Lin@meta.data[["sex"]] <- "F"
P08_Lin@meta.data[["age"]] <- "67"

saveRDS(P08_Lin, paste(dir,"/P08_Lin.rds",sep=""))

###MET01_Lin
MET01_Linfile <- "/blue/ferrallm/antonihalka/Pancreas/Lin/MET01/SRX8778055/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
MET01_Lin.counts <- Read10X(data.dir=MET01_Linfile)

MET01_Lin <- CreateSeuratObject(counts=MET01_Lin.counts, project="MET01_Lin")
MET01_Lin@meta.data[["stage"]] <- "4"
MET01_Lin@meta.data[["tissue type"]] <- "Metastasis"
MET01_Lin@meta.data[["site"]] <- "Liver"
MET01_Lin@meta.data[["treatment"]] <- "Unspecified"
MET01_Lin@meta.data[["treatment.type"]] <- "Unspecified"
MET01_Lin@meta.data[["sex"]] <- "M"
MET01_Lin@meta.data[["age"]] <- "43"

saveRDS(MET01_Lin, paste(dir,"/MET01_Lin.rds",sep=""))

###MET02_Lin
MET02_Linfile <- "/blue/ferrallm/antonihalka/Pancreas/Lin/MET02/SRX8778056/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
MET02_Lin.counts <- Read10X(data.dir=MET02_Linfile)

MET02_Lin <- CreateSeuratObject(counts=MET02_Lin.counts, project="MET02_Lin")
MET02_Lin@meta.data[["stage"]] <- "4"
MET02_Lin@meta.data[["tissue type"]] <- "Metastasis"
MET02_Lin@meta.data[["site"]] <- "Liver"
MET02_Lin@meta.data[["treatment"]] <- "Unspecified"
MET02_Lin@meta.data[["treatment.type"]] <- "Unspecified"
MET02_Lin@meta.data[["sex"]] <- "M"
MET02_Lin@meta.data[["age"]] <- "61"

saveRDS(MET02_Lin, paste(dir,"/MET02_Lin.rds",sep=""))

###MET04_Lin
MET04_Linfile <- "/blue/ferrallm/antonihalka/Pancreas/Lin/MET04/SRX8778058/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
MET04_Lin.counts <- Read10X(data.dir=MET04_Linfile)

MET04_Lin <- CreateSeuratObject(counts=MET04_Lin.counts, project="MET04_Lin")
MET04_Lin@meta.data[["stage"]] <- "4"
MET04_Lin@meta.data[["tissue type"]] <- "Metastasis"
MET04_Lin@meta.data[["site"]] <- "Liver"
MET04_Lin@meta.data[["treatment"]] <- "Unspecified"
MET04_Lin@meta.data[["treatment.type"]] <- "Unspecified"
MET04_Lin@meta.data[["sex"]] <- "F"
MET04_Lin@meta.data[["age"]] <- "47"

saveRDS(MET04_Lin, paste(dir,"/MET04_Lin.rds",sep=""))


###MET05_Lin
MET05_Linfile <- "/blue/ferrallm/antonihalka/Pancreas/Lin/MET05/SRX8778059/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
MET05_Lin.counts <- Read10X(data.dir=MET05_Linfile)

MET05_Lin <- CreateSeuratObject(counts=MET05_Lin.counts, project="MET05_Lin")
MET05_Lin@meta.data[["stage"]] <- "4"
MET05_Lin@meta.data[["tissue type"]] <- "Metastasis"
MET05_Lin@meta.data[["site"]] <- "Liver"
MET05_Lin@meta.data[["treatment"]] <- "Unspecified"
MET05_Lin@meta.data[["treatment.type"]] <- "Unspecified"
MET05_Lin@meta.data[["sex"]] <- "M"
MET05_Lin@meta.data[["age"]] <- "46"

saveRDS(MET05_Lin, paste(dir,"/MET05_Lin.rds",sep=""))


###MET06_Lin
MET06_Linfile <- "/blue/ferrallm/antonihalka/Pancreas/Lin/MET06/SRX8778060/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
MET06_Lin.counts <- Read10X(data.dir=MET06_Linfile)

MET06_Lin <- CreateSeuratObject(counts=MET06_Lin.counts, project="MET06_Lin")
MET06_Lin@meta.data[["stage"]] <- "4"
MET06_Lin@meta.data[["tissue type"]] <- "Metastasis"
MET06_Lin@meta.data[["site"]] <- "Liver"
MET06_Lin@meta.data[["treatment"]] <- "Unspecified"
MET06_Lin@meta.data[["treatment.type"]] <- "Unspecified"
MET06_Lin@meta.data[["sex"]] <- "M"
MET06_Lin@meta.data[["age"]] <- "67"

saveRDS(MET06_Lin, paste(dir,"/MET06_Lin.rds",sep=""))

###GSE148073#####################################

###Control1
Control1file <- "/blue/ferrallm/antonihalka/Pancreas/Fasolino_et_al/GSM4453619/Control1_HPAP022_scRNA-seq/SRX8053786/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Control1.counts <- Read10X(data.dir=Control1file)

# create a Seurat object based on the scRNA-seq data
Control1 <- CreateSeuratObject(counts=Control1.counts, project="Control1_GSE148073")
#add metadata
Control1@meta.data[["stage"]] <- "N/A"
Control1@meta.data[["tissue type"]] <- "Normal"
Control1@meta.data[["site"]] <- "Pancreas"
Control1@meta.data[["treatment"]] <- "N/A"
Control1@meta.data[["sex"]] <- "F"
Control1@meta.data[["age"]] <- "39"

#save rds file
saveRDS(Control1, paste(dir,"/Control1_GSE148073.rds",sep=""))

###Control3
Control3file <- "/blue/ferrallm/antonihalka/Pancreas/Fasolino_et_al/GSM4453621/Control3_HPAP035_scRNA-seq/SRX8053788/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Control3.counts <- Read10X(data.dir=Control3file)

# create a Seurat object based on the scRNA-seq data
Control3 <- CreateSeuratObject(counts=Control3.counts, project="Control3_GSE148073")
#add metadata
Control3@meta.data[["stage"]] <- "N/A"
Control3@meta.data[["tissue type"]] <- "Normal"
Control3@meta.data[["site"]] <- "Pancreas"
Control3@meta.data[["treatment"]] <- "N/A"
Control3@meta.data[["sex"]] <- "M"
Control3@meta.data[["age"]] <- "35"

#save rds file
saveRDS(Control3, paste(dir,"/Control3_GSE148073.rds",sep=""))

###Control5
Control5file <- "/blue/ferrallm/antonihalka/Pancreas/Fasolino_et_al/GSM4453623/Control5_HPAP037_scRNA-seq/SRX8053790/cellranger800_tutorial/outs/filtered_feature_bc_matrix"
Control5.counts <- Read10X(data.dir=Control5file)

# create a Seurat object based on the scRNA-seq data
Control5 <- CreateSeuratObject(counts=Control5.counts, project="Control5_GSE148073")
#add metadata
Control5@meta.data[["stage"]] <- "N/A"
Control5@meta.data[["tissue type"]] <- "Normal"
Control5@meta.data[["site"]] <- "Pancreas"
Control5@meta.data[["treatment"]] <- "N/A"
Control5@meta.data[["sex"]] <- "F"
Control5@meta.data[["age"]] <- "35"

#save rds file
saveRDS(Control5, paste(dir,"/Control5_GSE148073.rds",sep=""))

##### Merge all samples #####
pan <- merge(x = P01, y = list(
  P02,P03,P04,P05,P06,P07,P08,P09,P10,P11,P12,P13,P14,P15,P16,P17,P19,P20,P21,P24,P25,
  P03_Lin,P04_Lin,P05_Lin,P06_Lin,P07_Lin,P08_Lin,
  MET01_Lin,MET02_Lin,MET04_Lin,MET05_Lin,MET06_Lin,
  Control1,Control3,Control5),
  merge.data = TRUE,
  project = "pan_GSE205013_GSE154778_GSE148073")

# Join layers (raw counts)
pan <- JoinLayers(pan, assay = "RNA", new.layer = "counts")

# Save the merged object
saveRDS(pan, paste(dir, "/pancreas__GSE205013_GSE154778_GSE148073.rds",sep=""))


##################################################################
### QUALITY CONTROL - FILTER, NORMALIZE, SCALE DATA 
##################################################################
#Take a look at the data first
dim(pan)
# Example output: 36601 258723
summary(pan@meta.data$nFeature_RNA)

### FILTERING AND REMOVING MT- GENES 
## Determining how upper cutoff for nFeatures compares
nFeatUpper_pan <- mean(pan@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(pan@meta.data$nFeature_RNA, na.rm=TRUE)
nFeatLower_pan <- 450 # can filter out less cells later if it removes too many immune cells

# Make a violin plot to inspect nFeature_RNA distribution
# Extract metadata
meta <- pan@meta.data

# Make violin plot with ggplot2
pdf(paste(figdir, "/pancreas__GSE205013_GSE154778_GSE148073-ViolinPlot_nFeatureRNA.pdf", sep=""),
    width = 11, height = 6)

ggplot(meta, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  theme_classic() +
  labs(y = "nFeature_RNA", x = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# Store mitochondrial percentage in object meta data
pan <- PercentageFeatureSet(pan, pattern = "^MT-", col.name = "percent.mt")
perMitoUpper_pan <- 25 # upper threshold for mitochondrial content

# Update metadata after adding percent.mt
meta <- pan@meta.data
# Violin plot for percent.mt using ggplot2
pdf(paste(figdir, "/pancreas__GSE205013_GSE154778_GSE148073-ViolinPlot_PercentMito_ggplot2.pdf", sep=""),
    width = 11, height = 6)

ggplot(meta, aes(x = orig.ident, y = percent.mt)) +
  geom_violin(fill = "skyblue", color = "black") +
  geom_jitter(height = 0, width = 0.2, alpha = 0.2, size = 0.5) +
  theme_bw() +
  labs(title = "Percent Mitochondrial Genes by Sample",
       x = "Sample",
       y = "Percent.mt") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# Extract metadata
meta <- pan@meta.data

# Reshape data for ggplot2
library(reshape2)
meta_melt <- melt(meta, id.vars = "orig.ident", 
                  measure.vars = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  variable.name = "Metric", value.name = "Value")

# Create the violin plot
pdf(paste(figdir, "/pancreas__GSE205013_GSE154778_GSE148073-VlnPlot_nFeature+nCount+PercentMito_ggplot2.pdf", sep=""),
    width = 14, height = 6)

ggplot(meta_melt, aes(x = orig.ident, y = Value)) +
  geom_violin(fill = "skyblue", color = "black") +
  geom_jitter(height = 0, width = 0.2, alpha = 0.2, size = 0.5) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_bw() +
  labs(title = "QC Metrics by Sample", x = "Sample", y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# Scatter plots for QC
# Extract metadata
meta <- pan@meta.data

# Scatter plot: nCount_RNA vs percent.mt
pdf(paste(figdir, "/pancreas__GSE205013_GSE154778_GSE148073-FeatureScatter_nCountvMito_ggplot2.pdf", sep=""),
    width = 11, height = 6)

ggplot(meta, aes(x = nCount_RNA, y = percent.mt)) +
  geom_point(alpha = 0.3, size = 0.5, color = "darkblue") +
  theme_bw() +
  labs(title = "nCount_RNA vs Percent Mitochondrial Genes",
       x = "nCount_RNA",
       y = "Percent.mt")

dev.off()


# Scatter plot: nCount_RNA vs nFeature_RNA
pdf(paste(figdir, "/pancreas__GSE205013_GSE154778_GSE148073-FeatureScatter_nCountvnFeature_ggplot2.pdf", sep=""),
    width = 11, height = 6)

ggplot(meta, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(alpha = 0.3, size = 0.5, color = "darkgreen") +
  theme_bw() +
  labs(title = "nCount_RNA vs nFeature_RNA",
       x = "nCount_RNA",
       y = "nFeature_RNA")

dev.off()

## Original dataset dimensions
dim(pan)

## Filter based on features and mitochondrial content
pan <- subset(pan, subset = nFeature_RNA > nFeatLower_pan & nFeature_RNA < nFeatUpper_pan & percent.mt < perMitoUpper_pan)
saveRDS(pan, paste(dir, "/pancreas__GSE205013_GSE154778_GSE148073-PostFiltering_2026-02-22.rds", sep=""))
read()
pan <- readRDS("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/rds/pancreas__GSE205013_GSE154778_GSE148073-PostFiltering_2026-02-22.rds")
# Normalize with SCTransform
pan <- SCTransform(pan, verbose = FALSE)
dim(pan)

# Save normalized object
saveRDS(pan, paste(dir, "/pancreas__GSE205013_GSE154778_GSE148073-PostFiltering_2026-02-22.rds", sep=""))

#######################################################
### PCA, HARMONY, FINDNEIGHBORS, FINDCLUSTERS AND UMAP
#######################################################
### RUN PCA
#when running PCA after SCT, there's no need to specify VariableFeatures but it can be specified if needed 
pan <- RunPCA(pan, features=VariableFeatures(object=pan))

#PC_ 1 
#Positive:  COL3A1, COL1A2, COL1A1, SPP1, MGP, DCN, SFRP2, LUM, FTL, APOE 
#HLA-DRA, IGFBP7, SPARC, HBB, FN1, CCL4, CD74, B2M, C1QB, CALD1 
#HLA-DRB1, TAGLN, SPARCL1, COL6A3, APOC1, S100A9, C1QA, GNLY, VIM, IFI30 
#Negative:  PRSS2, PRSS1, REG1A, CELA3A, CPA1, CLPS, CTRB1, PLA2G1B, CPB1, CTRB2 
#CTRC, CELA3B, PNLIP, CPA2, SYCN, CELA2A, SPINK1, GP2, CEL, REG1B 
#AMY2A, PNLIPRP1, CELA2B, KLK1, PRSS3, INS, CTRL, REG3G, CUZD1, AMY2B 
#PC_ 2 
#Positive:  SPP1, FTL, APOE, HLA-DRA, CCL4, C1QB, CD74, APOC1, C1QA, HLA-DRB1 
#LYZ, IFI30, S100A9, GNLY, HLA-DPA1, C1QC, CCL5, TYROBP, HLA-DPB1, CCL3 
#SRGN, HLA-DQA1, NKG7, CXCL8, LCN2, TTR, HBB, FCER1G, CCL4L2, S100A8 
#Negative:  COL3A1, COL1A2, COL1A1, SFRP2, DCN, LUM, MGP, SPARC, IGFBP7, COL6A3 
#CALD1, TAGLN, FBLN1, FN1, CCDC80, IGFBP5, CTHRC1, C7, CCN2, SPARCL1 
#ASPN, POSTN, TIMP3, COL5A2, BGN, AEBP1, ACTA2, SFRP4, C11orf96, COL6A2 
#PC_ 3 
#Positive:  SPP1, APOE, FTL, C1QB, HLA-DRA, APOC1, C1QA, C1QC, IFI30, CD74 
#HLA-DRB1, HLA-DPA1, TYROBP, HLA-DPB1, HLA-DQA1, CTSB, CTSD, CCL18, S100A9, FTH1 
#GPNMB, FCER1G, CXCL8, AIF1, CD14, HMOX1, CCL3, CD68, IL1B, HLA-DQB1 
#Negative:  TTR, GCG, LCN2, INS, TFF2, TFF1, KRT19, HBB, PGC, TFF3 
#AGR2, S100A6, SCG5, TSPAN8, MUC6, CEACAM6, CHGB, PCSK1N, KRT18, TM4SF1 
#KRT8, MMP7, GNLY, CRYBA2, S100P, TM4SF4, OLFM4, MUC1, MUCL3, FXYD3 
#PC_ 4 
#Positive:  TTR, GCG, INS, SCG5, CHGB, CRYBA2, PCSK1N, CHGA, BEX1, FTL 
#FXYD2, SST, VGF, CLU, PTPRN, SCGN, SCG2, IAPP, PCSK2, CPE 
#SLC30A8, SCG3, GC, TM4SF4, MIR7-3HG, GPX3, NEUROD1, APOE, SPP1, STMN2 
#Negative:  HBB, GNLY, CCL4, REG1A, CCL5, HBA2, NKG7, TFF2, SPINK1, LCN2 
#TFF1, PRSS2, PRSS1, CD69, LYZ, HBA1, REG1B, GZMA, CXCR4, KRT19 
#PGC, PTPRC, GZMB, IL7R, S100A6, AGR2, S100A4, AMY2A, KLRB1, CPA1 
#PC_ 5 
#Positive:  TFF2, TFF1, LCN2, PGC, LYZ, REG1A, MUC6, S100A6, TFF3, KRT19 
#AGR2, SPP1, IFI27, SPINK1, TSPAN8, OLFM4, FTL, APOE, CEACAM6, MMP7 
#TM4SF1, KRT18, S100P, PRSS2, CLDN18, KRT8, MUC1, PIGR, MUC5AC, MUCL3 
#Negative:  CCL4, GNLY, HBB, CCL5, NKG7, CD69, GZMA, CXCR4, GZMB, HBA2 
#CCL4L2, PTPRC, KLRB1, IL7R, CST7, GZMK, KLRD1, DUSP2, PRF1, CD3D 
#CD2, SRGN, CTSW, TNFAIP3, HBA1, CCL3, CD52, BTG1, GZMH, CD7 

pdf(paste0(figdir, "/pancreas__GSE205013_GSE154778_GSE148073-PCA_DimPlot.pdf"), width = 11, height = 6)
DimPlot(pan, reduction = "pca", pt.size = 1)
dev.off()

# Determine dimensionality of dataset
pdf(paste0(figdir, "/pancreas__GSE205013_GSE154778_GSE148073-PCA_ElbowPlot.pdf"), width = 11, height = 6)
ElbowPlot(pan, ndims = 50) # from this plot, we can see that ~40 PCs capture most of the deviation
dev.off()

saveRDS(pan, paste0(dir, "/pancreas__GSE205013_GSE154778_GSE148073-PostPCA_2026-02-22.rds"))
pan <- readRDS("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/rds/pancreas__GSE205013_GSE154778_GSE148073-PostPCA_2026-02-22.rds")
### RUN HARMONY 
library('harmony')

panH <- RunHarmony(object=pan, group.by.vars="orig.ident")

pdf(paste0(figdir, "/pancreas__GSE205013_GSE154778_GSE148073-PCA_DimPlot_PostHarmony.pdf"), width = 11, height = 6)
DimPlot(panH, reduction = "harmony", pt.size = 1)
dev.off()

saveRDS(panH, paste0(dir, "/pancreas__GSE205013_GSE154778_GSE148073-PostHarmony_2026-02-23.rds"))

###FIND NEIGHBORS AND CLUSTERS
#for REGULAR dataset
pan <- FindNeighbors(pan, dims=1:50) #use all the dimensions, no need to reduce the dimensions for now
pan <- FindClusters(pan, resolution=seq(0.025,0.3, by=0.025)) # clusters

#Number of edges: 6231886

# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Maximum modularity in 10 random starts: 0.9690
#    Number of communities: 21
#    Elapsed time: 49 seconds
#    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#    
#    Number of nodes: 180740
#    Number of edges: 6231886
#    
#    Running Louvain algorithm...
#    0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|
#         Maximum modularity in 10 random starts: 0.9665
#       Number of communities: 26
#       Elapsed time: 49 seconds
#       Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#       
#       Number of nodes: 180740
#       Number of edges: 6231886
#       
#       Running Louvain algorithm...
#       0%   10   20   30   40   50   60   70   80   90   100%
#         [----|----|----|----|----|----|----|----|----|----|
#            **************************************************|
#            Maximum modularity in 10 random starts: 0.9642
#          Number of communities: 27
#          Elapsed time: 45 seconds
#          Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#          
#          Number of nodes: 180740
#          Number of edges: 6231886
#          
#          Running Louvain algorithm...
#          0%   10   20   30   40   50   60   70   80   90   100%
#            [----|----|----|----|----|----|----|----|----|----|
#               **************************************************|
#               Maximum modularity in 10 random starts: 0.9625
#             Number of communities: 26
#             Elapsed time: 50 seconds
#             Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#             
#             Number of nodes: 180740
#             Number of edges: 6231886
#             
#             Running Louvain algorithm...
#             0%   10   20   30   40   50   60   70   80   90   100%
#               [----|----|----|----|----|----|----|----|----|----|
#                  **************************************************|
#                  Maximum modularity in 10 random starts: 0.9609
#                Number of communities: 29
#                Elapsed time: 49 seconds
#                Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#                
#                Number of nodes: 180740
#                Number of edges: 6231886
#                
#                Running Louvain algorithm...
#                0%   10   20   30   40   50   60   70   80   90   100%
#                  [----|----|----|----|----|----|----|----|----|----|
#                     **************************************************|
#                     Maximum modularity in 10 random starts: 0.9596
#                   Number of communities: 32
#                   Elapsed time: 50 seconds

# Create the clustree plot to pick a resolution
pdf(paste0(figdir, "/pancreas__GSE205013_GSE154778_GSE148073-clustree.pdf"), width = 9, height = 9)
clustree(pan, prefix = "SCT_snn_res.")
dev.off() #0.075 or 0.1. O.075 looks better

resUse <- 0.075

pan$clusterResolution_0.075 <- as.factor(as.numeric(as.character(pan$SCT_snn_res.0.075)))

saveRDS(pan, paste(dir,"pancreas__GSE205013_GSE231525_GSE148073-PostNeighborANDClustering_2025-10-31",".rds", sep=""))


# For the HARMONY dataset
panH <- readRDS(paste0(dir, "pancreas__GSE205013_GSE154778_GSE148073-PostHarmony_2026-02-23.rds"))

# Compute neighbors using Harmony reduction
panH <- FindNeighbors(panH, dims = 1:50, reduction = "harmony")

# Run clustering across resolutions
panH <- FindClusters(panH, resolution = seq(0.025, 0.3, by = 0.025))


pdf(paste0(figdir, '/pancreas__GSE205013_GSE231525_GSE148073_withHarmony.pdf'),
    width = 9, height = 7)

clustree(panH, prefix = "SCT_snn_res.")

dev.off() #Both 0.075 and 0.1 look good

# Set resolutions to use for Harmony dataset
resUseH <- 0.075
panH$clusterResolution_0.075 <- as.factor(as.numeric(as.character(panH$SCT_snn_res.0.075)))

resUseH <- 0.1
panH$clusterResolution_0.1 <- as.factor(as.numeric(as.character(panH$SCT_snn_res.0.1)))

# Save updated Seurat object
saveRDS(panH,paste0(dir, "pancreas__GSE205013_GSE231525_GSE148073-PostNeighborANDClustering_withHarmony_2026-02-24.rds"))

### RUN UMAP
# For the REGULAR dataset
pan <- RunUMAP(pan, dims = 1:50) # UMAP based on PCA

# UMAP plots by cluster based on chosen resolution
pdf(paste0(figdir, "UMAP_res=", resUse, "-pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_clusters.pdf"), width = 15, height = 10)
DimPlot(pan, reduction = "umap", group.by = paste0("clusterResolution_", resUse))
dev.off()

# UMAP by sample
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_origIdent.pdf"), width = 15, height = 10)
DimPlot(pan, reduction = "umap", group.by = "orig.ident")
dev.off()

# UMAP by treatment
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_treatment.pdf"), width = 15, height = 10)
DimPlot(pan, reduction = "umap", group.by = "treatment")
dev.off()

# UMAP by tissue type
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_tissue_type.pdf"), width = 15, height = 10)
DimPlot(pan, reduction = "umap", group.by = "tissue.type")
dev.off()

# UMAP by site
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_site.pdf"), width = 15, height = 10)
DimPlot(pan, reduction = "umap", group.by = "site")
dev.off()

# UMAP by sex
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_sex.pdf"), width = 15, height = 10)
DimPlot(pan, reduction = "umap", group.by = "sex")
dev.off()

# UMAP by treatment.type
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_treatment_type.pdf"), width = 15, height = 10)
DimPlot(pan, reduction = "umap", group.by = "treatment.type")
dev.off()

# UMAP by stage
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_stage.pdf"), width = 15, height = 10)
DimPlot(pan, reduction = "umap", group.by = "stage")
dev.off()

# Save updated object
saveRDS(pan, paste0(dir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-PostUMAP_2026-02-24.rds"))

### For the HARMONY dataset
panH <- RunUMAP(panH, dims = 1:50, reduction = "harmony")

# UMAP plots by cluster for multiple resolutions
for(resUseH in c(0.075, 0.1)) {
  pdf(paste0(figdir, "UMAP_res=", resUseH, "-pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_clusters_withHarmony.pdf"), width = 15, height = 10)
  DimPlot(panH, reduction = "umap", group.by = paste0("clusterResolution_", resUseH))
  dev.off()
}

# UMAP by sample
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_origIdent_withHarmony.pdf"), width = 15, height = 10)
DimPlot(panH, reduction = "umap", group.by = "orig.ident")
dev.off()

# UMAP by treatment
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_treatment_withHarmony.pdf"), width = 15, height = 10)
DimPlot(panH, reduction = "umap", group.by = "treatment")
dev.off()

# UMAP by tissue type
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_tissue_type_withHarmony.pdf"), width = 15, height = 10)
DimPlot(panH, reduction = "umap", group.by = "tissue.type")
dev.off()

# UMAP by site
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_site_withHarmony.pdf"), width = 15, height = 10)
DimPlot(panH, reduction = "umap", group.by = "site")
dev.off()

# UMAP by sex
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_sex_withHarmony.pdf"), width = 15, height = 10)
DimPlot(panH, reduction = "umap", group.by = "sex")
dev.off()

# UMAP by treatement type
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_treatment_type_withHarmony.pdf"), width = 15, height = 10)
DimPlot(panH, reduction = "umap", group.by = "treatment.type")
dev.off()

# UMAP by stgae
pdf(paste0(figdir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-UMAP_stage_withHarmony.pdf"), width = 15, height = 10)
DimPlot(panH, reduction = "umap", group.by = "stage")
dev.off()

# Save updated Harmony object
saveRDS(panH, paste0(dir, "pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-PostUMAP_withHarmony_2026-02-24.rds"))




#################################################################################################################################
######################################Cell Type Assignment########################################################################

panH <- readRDS("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/rds/pancreas__GSE205013_GSE154778_GSE148073-JoinedLayers-PostUMAP_withHarmony_2026-02-24.rds")

#Check total # if cells
dim(panH)  #180740 cells

#Assign cells score based on genes to type them

#Acinar Cell genes
Acinar <- c('PRSS1', 'CTRB1','CTRB2','REG1B','CPA1','REG1A','SPINK1')
panH <- AddModuleScore(panH, features = list(Acinar), name = "Acinar Cells", search = TRUE)

#Activated Dendrtic Cells genes
Activated_DCs <- c("HLA-DRA","CCL22","LAMP3")
panH <- AddModuleScore(panH, features = list(Activated_DCs), name = "Activated Dendrtic Cells", search = TRUE)

# B cell genes
Bcell <- c('MS4A1','CD79A','CD79B','CD52', 'IGHM')
panH <- AddModuleScore(panH, features = list(Bcell), name = "B Cells", search = TRUE)

# Cycling Ductal genes
CyclingDuctal <- c('KRT19','ELF3','CENPF','MKI67')
panH <- AddModuleScore(panH, features = list(CyclingDuctal), name = "Cycling Ductal Cells", search = TRUE)

# Cycling Myeloid genes
CyclingMyeloid <- c('CENPF','MKI67','TYROBP','IFI30','STMN1','HMGB2','C1QB','C1QA')
panH <- AddModuleScore(panH, features = list(CyclingMyeloid), name = "Cycling Myeloid", search = TRUE)

# Cycling TNK genes (this detonates T and NK cells)
CyclingTNK <- c('CENPF','MKI67','CCL5','NKG7','STMN1','HMGB2')
panH <- AddModuleScore(panH, features = list(CyclingTNK), name = "Cycling TNK", search = TRUE)

# Ductal Cell 1 genes
Ductal1 <- c('MBP','CFTR','MMP7','KRT18','KRT8','TFF1','KRT19')
panH <- AddModuleScore(panH, features = list(Ductal1), name = "Ductal Cell 1", search = TRUE)

# Ductal Cell 2 genes
Ductal2 <- c('KRT19','KRT7','TSPAN8','SLPI','MMP7','SPP1','CLU')
panH <- AddModuleScore(panH, features = list(Ductal2), name = "Ductal Cell 2", search = TRUE)

# Endocrine cell genes
Endocrine <- c('CHGB','CHGA','INS','IAPP','TTR','ERO1B')
panH <- AddModuleScore(panH, features = list(Endocrine), name = "Endocrine Cells", search = TRUE)

# Endothelial cell genes
Endothelial <- c('CDH5','PLVAP','VWF','CLDN5','RAMP2')
panH <- AddModuleScore(panH, features = list(Endothelial), name = "Endothelial Cells", search = TRUE)

# Fibroblast genes
Fibroblast <- c('LUM','DCN','COL1A1','COL1A2','TAGLN','ACTA2','COL3A1')
panH <- AddModuleScore(panH, features = list(Fibroblast), name = "Fibroblasts", search = TRUE)

# Macrophage genes
Macrophage <- c('AIF1','FCGR1A','CD14','CD68') #Note that FCGR1A is the name for CD64
panH <- AddModuleScore(panH, features = list(Macrophage), name = "Macrophages", search = TRUE)

# Mast cell genes
Mast <- c('TPSAB1','TPSB2','CPA3')
panH <- AddModuleScore(panH, features = list(Mast), name = "Mast Cells", search = TRUE)

# Myeloid genes
Myeloid <- c('TYROBP','IFI30','LYZ','APOE','HLA-DRA','C1QA','CD14')
panH <- AddModuleScore(panH, features = list(Myeloid), name = "Myeloid Cells", search = TRUE)

# pDC genes
pDC <- c('IRF7','GZMB')
panH <- AddModuleScore(panH, features = list(pDC), name = "pDCs", search = TRUE)

# Pericyte genes
Pericytes <- c('TAGLN','ACTA2')
panH <- AddModuleScore(panH, features = list(Pericytes), name = "Pericytes", search = TRUE)

# Perivascular cell genes
Perivascular <- c('IGFBP7','ACTA2','RGS5')
panH <- AddModuleScore(panH, features = list(Perivascular), name = "Perivascular Cells", search = TRUE)

# Plasma cell genes
Plasma <- c('IGKC','IGHG1','IGJ','IGLL5','CD79A') #Could not find IGLL5
panH <- AddModuleScore(panH, features = list(Plasma), name = "Plasma Cells", search = TRUE)

# Red Blood Cell genes
RBC <- c('HBA2','HBB','MMP7')
panH <- AddModuleScore(panH, features = list(RBC), name = "Red Blood Cells", search = TRUE)

# Stellate cell genes
Stellate <- c('RGS5','ACTA2','PDGFRB','ADIRF')
panH <- AddModuleScore(panH, features = list(Stellate), name = "Stellate Cells", search = TRUE)

# T Cell genes
T_Cells <- c('CD3D','CD3E','CD4','CD8')
panH <- AddModuleScore(panH, features = list(T_Cells), name = "T Cells", search = TRUE)

# NK Cell genes
NK_Cells <- c('NKG7','CD2','CCL5','NCAM1')  # NCAM1 = CD56
panH <- AddModuleScore(panH, features = list(NK_Cells), name = "NK Cells", search = TRUE)

saveRDS(panH, paste(dir,"/pancreas__GSE205013_GSE154778_GSE148073_withHarmony_PostCellTyping_2026-02-25", ".rds", sep=""))

##Group Cells
cell_type_cols <- c('Acinar Cells1', 'Activated Dendrtic Cells1', 'B Cells1', 'Cycling Ductal Cells1', 'Cycling TNK1', 'Ductal Cell 11', 'Ductal Cell 21',
                    'Endocrine Cells1', 'Endothelial Cells1', 'Fibroblasts1','Myeloid Cells1', 'Cycling Myeloid1', 'Mast Cells1', 'Macrophages1',
                    'NK Cells1','pDCs1', 'Pericytes1', 'Perivascular Cells1', 'Plasma Cells1','Red Blood Cells1', 'T Cells1', 'Stellate Cells1')

panH$AssignedCell <- apply(panH@meta.data[, cell_type_cols, drop = FALSE], 1, function(x) cell_type_cols[which.max(x)])
panH$AssignedCell <- gsub("1$", "", panH$AssignedCell)  

table(panH$AssignedCell)

#Grouping
panH$CellGroup <- panH$AssignedCell

#Ductal
panH$CellGroup[panH$AssignedCell %in% c("Ductal Cell 1", "Ductal Cell 2", "Cycling Ductal Cells")] <- "Ductal"

#Myeloid
panH$CellGroup[panH$AssignedCell %in% c("Myeloid Cells", "Cycling Myeloid", "Macrophages", "Activated Dendritic Cells")] <- "Myeloid"

#T/NK
panH$CellGroup[panH$AssignedCell %in% c("NK Cells", "Cycling TNK", "T Cells")] <- "T/NK"

#Perivascular
panH$CellGroup[panH$AssignedCell %in% c("Pericytes", "Perivascular Cells")] <- "Perivascular"


# Fix typo in tissue.type
panH$tissue.type[panH$tissue.type == "Metastatsis"] <- "Metastasis"

group_table <- table(panH$CellGroup, panH$tissue.type)
head(group_table)
write.csv(group_table, "/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/tables/cellgroup_counts.csv", row.names = TRUE)

saveRDS(panH, paste(dir,"/pancreas__GSE205013_GSE154778_GSE148073_withHarmony_CellGroupingsTyping_2026-02-25", ".rds", sep=""))

########################
###Regulated Pathways###
########################

panH <- readRDS("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/rds/pancreas__GSE205013_GSE154778_GSE148073_withHarmony_CellGroupingsTyping_2026-02-25.rds")
Idents(panH) <- "tissue.type"

# Run FindMarkers
prim_v_norm.markers <- FindMarkers(panH, ident.1 = "Primary", ident.2 = "Normal")
prim_v_norm.markers$gene <- rownames(prim_v_norm.markers)
sig_prim_norm <- prim_v_norm.markers %>% filter(p_val_adj <= 0.05)

meta_v_norm.markers <- FindMarkers(panH, ident.1 = "Metastasis", ident.2 = "Normal")
meta_v_norm.markers$gene <- rownames(meta_v_norm.markers)
sig_meta_norm <- meta_v_norm.markers %>% filter(p_val_adj <= 0.05)

prim_v_meta.markers <- FindMarkers(panH, ident.1 = "Primary", ident.2 = "Metastasis")
prim_v_meta.markers$gene <- rownames(prim_v_meta.markers)
sig_prim_meta <- prim_v_meta.markers %>% filter(p_val_adj <= 0.05)

### Upregulated Pathways
# Primary vs Normal
sig_prim_norm_up <- sig_prim_norm %>% filter(avg_log2FC > 1)
sig_prim_norm_up_map <- bitr(sig_prim_norm_up$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
bg_sig_prim_norm_map <- bitr(rownames(panH), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

m_Hpath <- msigdbr(species="Homo sapiens", category="H") %>% dplyr::select(gs_name, entrez_gene)
msig_Hpath_pvn_upgenes <- enricher(sig_prim_norm_up_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_prim_norm_map$ENTREZID)

msig_Hpath_pvn_upgenes_df <- msig_Hpath_pvn_upgenes@result
write.csv(msig_Hpath_pvn_upgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Prim_v_Norm.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Prim_v_Norm.pdf'), width=12, height=12)
barplot(msig_Hpath_pvn_upgenes, showCategory=20)
dev.off()

# Metastasis vs Normal
sig_meta_norm_up <- sig_meta_norm %>% filter(avg_log2FC > 1)
sig_meta_norm_up_map <- bitr(sig_meta_norm_up$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
bg_sig_meta_norm_map <- bitr(rownames(panH), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

msig_Hpath_mvn_upgenes <- enricher(sig_meta_norm_up_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_meta_norm_map$ENTREZID)
msig_Hpath_mvn_upgenes_df <- msig_Hpath_mvn_upgenes@result
write.csv(msig_Hpath_mvn_upgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Meta_v_Norm.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Meta_v_Norm.pdf'), width=12, height=12)
barplot(msig_Hpath_mvn_upgenes, showCategory=20)
dev.off()

# Primary vs Metastasis
sig_prim_meta_up <- sig_prim_meta %>% filter(avg_log2FC > 1)
sig_prim_meta_up_map <- bitr(sig_prim_meta_up$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
bg_sig_prim_meta_map <- bitr(rownames(panH), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

msig_Hpath_pvm_upgenes <- enricher(sig_prim_meta_up_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_prim_meta_map$ENTREZID)
msig_Hpath_pvm_upgenes_df <- msig_Hpath_pvm_upgenes@result
write.csv(msig_Hpath_pvm_upgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Prim_v_Meta.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Prim_v_Meta.pdf'), width=12, height=12)
barplot(msig_Hpath_pvm_upgenes, showCategory=20)
dev.off()

### Downregulated Pathways
# Primary vs Normal
sig_prim_norm_down <- sig_prim_norm %>% filter(avg_log2FC <= -1)
sig_prim_norm_down_map <- bitr(sig_prim_norm_down$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
msig_Hpath_pvn_downgenes <- enricher(sig_prim_norm_down_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_prim_norm_map$ENTREZID)
msig_Hpath_pvn_downgenes_df <- msig_Hpath_pvn_downgenes@result
write.csv(msig_Hpath_pvn_downgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Prim_v_Norm.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Prim_v_Norm.pdf'), width=12, height=12)
barplot(msig_Hpath_pvn_downgenes, showCategory=20)
dev.off()

# Metastasis vs Normal
sig_meta_norm_down <- sig_meta_norm %>% filter(avg_log2FC <= -1)
sig_meta_norm_down_map <- bitr(sig_meta_norm_down$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
msig_Hpath_mvn_downgenes <- enricher(sig_meta_norm_down_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_meta_norm_map$ENTREZID)
msig_Hpath_mvn_downgenes_df <- msig_Hpath_mvn_downgenes@result
write.csv(msig_Hpath_mvn_downgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Meta_v_Norm.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Meta_v_Norm.pdf'), width=12, height=12)
barplot(msig_Hpath_mvn_downgenes, showCategory=20)
dev.off()

# Primary vs Metastasis
sig_prim_meta_down <- sig_prim_meta %>% filter(avg_log2FC <= -1)
sig_prim_meta_down_map <- bitr(sig_prim_meta_down$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
msig_Hpath_pvm_downgenes <- enricher(sig_prim_meta_down_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_prim_meta_map$ENTREZID)
msig_Hpath_pvm_downgenes_df <- msig_Hpath_pvm_downgenes@result
write.csv(msig_Hpath_pvm_downgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Prim_v_Meta.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Prim_v_Meta.pdf'), width=12, height=12)
barplot(msig_Hpath_pvm_downgenes, showCategory=20)
dev.off()

# Extract Significant Genes
extract_pathway_genes <- function(enrich_obj, comparison_label) {
  if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) {
    return(data.frame(Pathway=character(), adj_pval=numeric(), GeneCount=integer(),
                      ENTREZID=character(), SYMBOL=character(), Comparison=character(),
                      stringsAsFactors=FALSE))
  }
  df <- as_tibble(enrich_obj@result) %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::select(Pathway=Description, adj_pval=p.adjust, GeneCount=Count, ENTREZID=geneID) %>%
    tidyr::separate_rows(ENTREZID, sep="/") %>%
    dplyr::mutate(Comparison=comparison_label)
  
  gene_map <- bitr(df$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb=org.Hs.eg.db)
  df <- df %>% dplyr::left_join(gene_map, by=c("ENTREZID"="ENTREZID"))
  return(df)
}

# Apply extraction
pvn_up_genes   <- extract_pathway_genes(msig_Hpath_pvn_upgenes, "PVN_Up")
pvn_down_genes <- extract_pathway_genes(msig_Hpath_pvn_downgenes, "PVN_Down")
mvn_up_genes   <- extract_pathway_genes(msig_Hpath_mvn_upgenes, "MVN_Up")
mvn_down_genes <- extract_pathway_genes(msig_Hpath_mvn_downgenes, "MVN_Down")
pvm_up_genes   <- extract_pathway_genes(msig_Hpath_pvm_upgenes, "PVM_Up")
pvm_down_genes <- extract_pathway_genes(msig_Hpath_pvm_downgenes, "PVM_Down")

all_pathway_genes <- dplyr::bind_rows(
  pvn_up_genes, pvn_down_genes,
  mvn_up_genes, mvn_down_genes,
  pvm_up_genes, pvm_down_genes
)

write.csv(all_pathway_genes,
          paste0(figdir, "Hallmark_Pathway_Driver_Genes_All_Comparisons_Spring2026.csv"),
          row.names=FALSE)

# Extract Significant Downregulated Genes
pvn_down_genes   <- extract_pathway_genes(msig_Hpath_pvn_downgenes, "PVN_Down")
mvn_down_genes   <- extract_pathway_genes(msig_Hpath_mvn_downgenes, "MVN_Down")
pvm_down_genes   <- extract_pathway_genes(msig_Hpath_pvm_downgenes, "PVM_Down")

all_downregulated_pathway_genes <- dplyr::bind_rows(
  pvn_down_genes,
  mvn_down_genes,
  pvm_down_genes
)

write.csv(all_downregulated_pathway_genes,
          paste0(figdir, "Hallmark_Pathway_Driver_Genes_All_Downregulated_Spring2026.csv"),
          row.names=FALSE)

##################################
##################################
######By Disease Progression######
##################################

panH <- readRDS("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/rds/pancreas__GSE205013_GSE154778_GSE148073_withHarmony_CellGroupingsTyping_2026-02-25.rds")

# Confirm stage values first
unique(panH$stage)

# Initialize
panH$progression_group <- NA_character_

# Assign groups explicitly
panH$progression_group[panH$stage %in% c(1, 2)] <- "Early"
panH$progression_group[panH$stage %in% c(3, 4)] <- "Advanced"
panH$progression_group[panH$stage %in% c('N/A')] <- "Control"

saveRDS(panH, paste(dir,"/pancreas__GSE205013_GSE154778_GSE148073_withHarmony_CellGroupingsTyping_VS_Progression_2026-03-02", ".rds", sep=""))

#####Cell groupings
group_table <- table(panH$CellGroup, panH$progression_group)
head(group_table)
write.csv(group_table, "/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/tables/cellgroup_vs_progression_counts.csv", row.names = TRUE)


########################
### Regulated Pathways #
########################
panH <- readRDS("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/rds/pancreas__GSE205013_GSE154778_GSE148073_withHarmony_CellGroupingsTyping_VS_Progression_2026-03-02.rds")

Idents(panH) <- "progression_group"

# Run FindMarkers
early_v_control.markers <- FindMarkers(panH, ident.1 = "Early", ident.2 = "Control")
early_v_control.markers$gene <- rownames(early_v_control.markers)
sig_early_control <- early_v_control.markers %>% filter(p_val_adj <= 0.05)

adv_v_control.markers <- FindMarkers(panH, ident.1 = "Advanced", ident.2 = "Control")
adv_v_control.markers$gene <- rownames(adv_v_control.markers)
sig_adv_control <- adv_v_control.markers %>% filter(p_val_adj <= 0.05)

early_v_adv.markers <- FindMarkers(panH, ident.1 = "Early", ident.2 = "Advanced")
early_v_adv.markers$gene <- rownames(early_v_adv.markers)
sig_early_adv <- early_v_adv.markers %>% filter(p_val_adj <= 0.05)

### Upregulated Pathways

# Early vs Control
sig_early_control_up <- sig_early_control %>% filter(avg_log2FC > 1)
sig_early_control_up_map <- bitr(sig_early_control_up$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
bg_sig_early_control_map <- bitr(rownames(panH), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

m_Hpath <- msigdbr(species="Homo sapiens", category="H") %>% dplyr::select(gs_name, entrez_gene)
msig_Hpath_evc_upgenes <- enricher(sig_early_control_up_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_early_control_map$ENTREZID)

msig_Hpath_evc_upgenes_df <- msig_Hpath_evc_upgenes@result
write.csv(msig_Hpath_evc_upgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Early_v_Control.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Early_v_Control.pdf'), width=12, height=12)
barplot(msig_Hpath_evc_upgenes, showCategory=20)
dev.off()

# Advanced vs Control
sig_adv_control_up <- sig_adv_control %>% filter(avg_log2FC > 1)
sig_adv_control_up_map <- bitr(sig_adv_control_up$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
bg_sig_adv_control_map <- bitr(rownames(panH), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

msig_Hpath_avc_upgenes <- enricher(sig_adv_control_up_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_adv_control_map$ENTREZID)
msig_Hpath_avc_upgenes_df <- msig_Hpath_avc_upgenes@result
write.csv(msig_Hpath_avc_upgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Advanced_v_Control.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Advanced_v_Control.pdf'), width=12, height=12)
barplot(msig_Hpath_avc_upgenes, showCategory=20)
dev.off()

# Early vs Advanced
sig_early_adv_up <- sig_early_adv %>% filter(avg_log2FC > 1)
sig_early_adv_up_map <- bitr(sig_early_adv_up$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
bg_sig_early_adv_map <- bitr(rownames(panH), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

msig_Hpath_eva_upgenes <- enricher(sig_early_adv_up_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_early_adv_map$ENTREZID)
msig_Hpath_eva_upgenes_df <- msig_Hpath_eva_upgenes@result
write.csv(msig_Hpath_eva_upgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Early_v_Advanced.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Early_v_Advanced.pdf'), width=12, height=12)
barplot(msig_Hpath_eva_upgenes, showCategory=20)
dev.off()

### Downregulated Pathways

# Early vs Control
sig_early_control_down <- sig_early_control %>% filter(avg_log2FC <= -1)
sig_early_control_down_map <- bitr(sig_early_control_down$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
msig_Hpath_evc_downgenes <- enricher(sig_early_control_down_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_early_control_map$ENTREZID)
msig_Hpath_evc_downgenes_df <- msig_Hpath_evc_downgenes@result
write.csv(msig_Hpath_evc_downgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Early_v_Control.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Early_v_Control.pdf'), width=12, height=12)
barplot(msig_Hpath_evc_downgenes, showCategory=20)
dev.off()

# Advanced vs Control
sig_adv_control_down <- sig_adv_control %>% filter(avg_log2FC <= -1)
sig_adv_control_down_map <- bitr(sig_adv_control_down$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
msig_Hpath_avc_downgenes <- enricher(sig_adv_control_down_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_adv_control_map$ENTREZID)
msig_Hpath_avc_downgenes_df <- msig_Hpath_avc_downgenes@result
write.csv(msig_Hpath_avc_downgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Advanced_v_Control.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Advanced_v_Control.pdf'), width=12, height=12)
barplot(msig_Hpath_avc_downgenes, showCategory=20)
dev.off()

# Early vs Advanced
sig_early_adv_down <- sig_early_adv %>% filter(avg_log2FC <= -1)
sig_early_adv_down_map <- bitr(sig_early_adv_down$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
msig_Hpath_eva_downgenes <- enricher(sig_early_adv_down_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_early_adv_map$ENTREZID)
msig_Hpath_eva_downgenes_df <- msig_Hpath_eva_downgenes@result
write.csv(msig_Hpath_eva_downgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Early_v_Advanced.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Early_v_Advanced.pdf'), width=12, height=12)
barplot(msig_Hpath_eva_downgenes, showCategory=20)
dev.off()

# Extract downregulated genes from each progression comparison
evc_down_genes <- extract_pathway_genes(msig_Hpath_evc_downgenes, "Early_v_Control_Down")
avc_down_genes <- extract_pathway_genes(msig_Hpath_avc_downgenes, "Advanced_v_Control_Down")
eva_down_genes <- extract_pathway_genes(msig_Hpath_eva_downgenes, "Early_v_Advanced_Down")

# Combine all downregulated genes
all_downregulated_pathway_genes <- dplyr::bind_rows(
  evc_down_genes,
  avc_down_genes,
  eva_down_genes
)

# Save combined CSV
write.csv(all_downregulated_pathway_genes,
          paste0(figdir, "mSigHpathways_All_Downregulated_Pathway_Genes_By_Progression_Spring2026.csv"),
          row.names = FALSE)

##################################
######By Treatment################
##################################

panH <- readRDS("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/rds/pancreas__GSE205013_GSE154778_GSE148073_withHarmony_CellGroupingsTyping_VS_Progression_2026-03-02.rds")

#####Cell groupings
group_table <- table(panH$CellGroup, panH$treatment) ##manually remove or ignore "Unspecified" column from output table
head(group_table)
write.csv(group_table, "/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/tables/cellgroup_vs_treatment_counts.csv", row.names = TRUE)

#regulated pathways

Idents(panH) <- "treatment"

# Run FindMarkers
treated_v_untreated.markers <- FindMarkers(panH, ident.1 = "Treated", ident.2 = "Untreated")
treated_v_untreated.markers$gene <- rownames(treated_v_untreated.markers)
sig_treated_untreated <- treated_v_untreated.markers %>% filter(p_val_adj <= 0.05)

treated_v_control.markers <- FindMarkers(panH, ident.1 = "Treated", ident.2 = "N/A")
treated_v_control.markers$gene <- rownames(treated_v_control.markers)
sig_treated_v_control <- treated_v_control.markers %>% filter(p_val_adj <= 0.05)

untreated_v_control.markers <- FindMarkers(panH, ident.1 = "Untreated", ident.2 = "N/A")
untreated_v_control.markers$gene <- rownames(untreated_v_control.markers)
sig_untreated_v_control <- untreated_v_control.markers %>% filter(p_val_adj <= 0.05)

### Upregulated Pathways

# Treated vs Untreated
sig_treated_untreated_up <- sig_treated_untreated %>% filter(avg_log2FC > 1)
sig_treated_untreated_up_map <- bitr(sig_treated_untreated_up$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
bg_sig_treated_untreated_map <- bitr(rownames(panH), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

m_Hpath <- msigdbr(species="Homo sapiens", category="H") %>% dplyr::select(gs_name, entrez_gene)
msig_Hpath_tvu_upgenes <- enricher(sig_treated_untreated_up_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_treated_untreated_map$ENTREZID)

msig_Hpath_tvu_upgenes_df <- msig_Hpath_tvu_upgenes@result
write.csv(msig_Hpath_tvu_upgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Treated_v_Untreated.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Treated_v_Untreated.pdf'), width=12, height=12)
barplot(msig_Hpath_tvu_upgenes, showCategory=20)
dev.off()

# Treated vs Control
sig_treated_v_control_up <- sig_treated_v_control %>% filter(avg_log2FC > 1)
sig_treated_v_control_up_map <- bitr(sig_treated_v_control_up$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
bg_sig_treated_v_control_map <- bitr(rownames(panH), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

msig_Hpath_tvc_upgenes <- enricher(sig_treated_v_control_up_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_treated_v_control_map$ENTREZID)
msig_Hpath_tvc_upgenes_df <- msig_Hpath_tvc_upgenes@result
write.csv(msig_Hpath_tvc_upgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Treated_v_Control.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Treated_v_Control.pdf'), width=12, height=12)
barplot(msig_Hpath_tvc_upgenes, showCategory=20)
dev.off()

# Untreated vs Control
sig_untreated_v_control_up <- sig_untreated_v_control %>% filter(avg_log2FC > 1)
sig_untreated_v_control_up_map <- bitr(sig_untreated_v_control_up$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
bg_sig_untreated_v_control_map <- bitr(rownames(panH), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

msig_Hpath_uvc_upgenes <- enricher(sig_untreated_v_control_up_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_untreated_v_control_map$ENTREZID)
msig_Hpath_uvc_upgenes_df <- msig_Hpath_uvc_upgenes@result
write.csv(msig_Hpath_uvc_upgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Untreated_v_Control.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Untreated_v_Control.pdf'), width=12, height=12)
barplot(msig_Hpath_uvc_upgenes, showCategory=20)
dev.off()

### Downregulated Pathways

# Treated vs Untreated
sig_treated_untreated_down <- sig_treated_untreated %>% filter(avg_log2FC <= -1)
sig_treated_untreated_down_map <- bitr(sig_treated_untreated_down$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
msig_Hpath_tvu_downgenes <- enricher(sig_treated_untreated_down_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_treated_untreated_map$ENTREZID)
msig_Hpath_tvu_downgenes_df <- msig_Hpath_tvu_downgenes@result
write.csv(msig_Hpath_tvu_downgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Treated_v_Untreated.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Treated_v_Untreated.pdf'), width=12, height=12)
barplot(msig_Hpath_tvu_downgenes, showCategory=20)
dev.off()

# Treated vs Control
sig_treated_v_control_down <- sig_treated_v_control %>% filter(avg_log2FC <= -1)
sig_treated_v_control_down_map <- bitr(sig_treated_v_control_down$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
msig_Hpath_tvc_downgenes <- enricher(sig_treated_v_control_down_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_treated_v_control_map$ENTREZID)
msig_Hpath_tvc_downgenes_df <- msig_Hpath_tvc_downgenes@result
write.csv(msig_Hpath_tvc_downgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Treated_v_Control.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Treated_v_Control.pdf'), width=12, height=12)
barplot(msig_Hpath_tvc_downgenes, showCategory=20)
dev.off()

# Untreated vs Control
sig_untreated_v_control_down <- sig_untreated_v_control %>% filter(avg_log2FC <= -1)
sig_untreated_v_control_down_map <- bitr(sig_untreated_v_control_down$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
msig_Hpath_uvc_downgenes <- enricher(sig_untreated_v_control_down_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_untreated_v_control_map$ENTREZID)
msig_Hpath_uvc_downgenes_df <- msig_Hpath_uvc_downgenes@result
write.csv(msig_Hpath_uvc_downgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Untreated_v_Control.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Untreated_v_Control.pdf'), width=12, height=12)
barplot(msig_Hpath_uvc_downgenes, showCategory=20)
dev.off()

# Untreated vs Control
sig_untreated_v_control_down <- sig_untreated_v_control %>% filter(avg_log2FC <= -1)
sig_untreated_v_control_down_map <- bitr(sig_untreated_v_control_down$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
msig_Hpath_uvc_downgenes <- enricher(sig_untreated_v_control_down_map$ENTREZID, TERM2GENE=m_Hpath, universe=bg_sig_untreated_v_control_map$ENTREZID)
msig_Hpath_uvc_downgenes_df <- msig_Hpath_uvc_downgenes@result
write.csv(msig_Hpath_uvc_downgenes_df, paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Untreated_v_Control.csv'))
pdf(paste0(figdir, 'mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Untreated_v_Control.pdf'), width=12, height=12)
barplot(msig_Hpath_uvc_downgenes, showCategory=20)
dev.off()

# Extract downregulated genes from each treatment comparison
tvu_down_genes <- extract_pathway_genes(msig_Hpath_tvu_downgenes, "Treated_v_Untreated_Down")
tvc_down_genes <- extract_pathway_genes(msig_Hpath_tvc_downgenes, "Treated_v_Control_Down")
uvc_down_genes <- extract_pathway_genes(msig_Hpath_uvc_downgenes, "Untreated_v_Control_Down")

# Combine all downregulated genes into a single table
all_downregulated_treatment_genes <- dplyr::bind_rows(
  tvu_down_genes,
  tvc_down_genes,
  uvc_down_genes
)

# Save combined CSV
write.csv(all_downregulated_treatment_genes,
          paste0(figdir, "mSigHpathways_All_Downregulated_Pathway_Genes_By_Treatment_Spring2026.csv"),
          row.names = FALSE)
saveRDS(panH, paste(dir,"/pancreas__GSE205013_GSE154778_GSE148073_withHarmony_2026-03-26", ".rds", sep=""))

panH <- readRDS("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/rds/pancreas__GSE205013_GSE154778_GSE148073_withHarmony_2026-03-26.rds")
#########################################################
#### extract data for COMET in metastatic vs primary ####
#########################################################

figdir2 <- "/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/COMET/"

Idents(panH) <- "sample_type"

# DEG between metastatic and primary cells
met_vs_primary.markers <- FindMarkers(panH, ident.1 = "Metastatic", ident.2 = "Primary")
met_vs_primary.markers$gene <- rownames(met_vs_primary.markers)

# significant genes
sig_met_primary_genes <- met_vs_primary.markers %>% 
  filter(p_val_adj <= 0.05)

# upregulated genes in metastatic
sig_met_primary_up <- sig_met_primary_genes %>% 
  filter(avg_log2FC > 1)

# normalized expression matrix
panH_counts_matrix <- GetAssayData(panH, slot = "data", assay = "SCT")

# metastatic cells
met_cells <- panH$sample_type == "Metastatic"
met.data <- panH_counts_matrix[,met_cells]

# UMAP
panH_embeddings <- Embeddings(panH, reduction = "umap")
met.umap <- panH_embeddings[met_cells,]

# primary cells
primary_cells <- panH$sample_type == "Primary"
primary.data <- panH_counts_matrix[,primary_cells]
primary.umap <- panH_embeddings[primary_cells,]

# COMET max = 65k cells
nSamples <- 64999 - dim(met.data)[2]

# cluster labels
met.cluster <- data.frame("X" = rep(1, dim(met.data)[2]))
primary.cluster <- data.frame("X" = rep(0, nSamples))

# convert to dataframe
met.data.df <- data.frame(met.data[sig_met_primary_up$gene,])
primary.data.df <- data.frame(primary.data[sig_met_primary_up$gene,])

# downsample primary cells
cellsKeep <- sample(length(primary.data.df), nSamples)

primary.data.DS <- primary.data.df[,cellsKeep]
primary.umap.DS <- data.frame(primary.umap[cellsKeep,])

# combine
markers.df <- cbind(met.data.df, primary.data.DS)
umap.df <- rbind(met.umap, primary.umap.DS)
clusters.df <- rbind(met.cluster, primary.cluster)

row.names(clusters.df) <- row.names(umap.df)
colnames(markers.df) <- row.names(umap.df)

# write COMET input files
write.table(markers.df, paste0(figdir2,"met_vs_primary_markers.txt"),
            sep="\t",row.names=TRUE, quote=FALSE)

write.table(umap.df, paste0(figdir2,"met_vs_primary_umap.txt"),
            sep="\t",row.names=TRUE,col.names=FALSE, quote=FALSE)

write.table(clusters.df, paste0(figdir2,"met_vs_primary_clusters.txt"),
            sep="\t",row.names=TRUE,col.names=FALSE, quote=FALSE)

write.table(sig_met_primary_up$gene,
            paste0(figdir2,"met_vs_primary_genes.txt"),
            sep="\t",row.names=FALSE,col.names=FALSE, quote=FALSE)

# checks
intersect(colnames(markers.df), rownames(clusters.df))

genes <- intersect(sig_met_primary_up$gene, rownames(markers.df))
all <- intersect(rownames(clusters.df), rownames(umap.df))


#########################################################
#### extract data for COMET in treated vs untreated ####
#########################################################

figdir3 <- "/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/COMET_TvU/"

Idents(panH) <- "treatment"

# DEG between treated and untreated cells
treated_vs_untreated.markers <- FindMarkers(panH, ident.1 = "Treated", ident.2 = "Untreated")
treated_vs_untreated.markers$gene <- rownames(treated_vs_untreated.markers)

# significant genes
sig_treated_untreated_genes <- treated_vs_untreated.markers %>% 
  filter(p_val_adj <= 0.05)

# upregulated genes in treated
sig_treated_untreated_up <- sig_treated_untreated_genes %>% 
  filter(avg_log2FC > 1)

# normalized expression matrix
panH_counts_matrix <- GetAssayData(panH, layer = "data", assay = "SCT")

# treated cells
treated_cells <- panH$treatment == "Treated"
treated.data <- panH_counts_matrix[,treated_cells]

# UMAP
panH_embeddings <- Embeddings(panH, reduction = "umap")
treated.umap <- panH_embeddings[treated_cells,]

# untreated cells
untreated_cells <- panH$treatment == "Untreated"
untreated.data <- panH_counts_matrix[,untreated_cells]
untreated.umap <- panH_embeddings[untreated_cells,]

# COMET max = 65k cells, balanced between groups
set.seed(1)  # for reproducibility
cells_per_group <- min(32500, ncol(treated.data), ncol(untreated.data))

# convert to dataframe
treated.data.df <- data.frame(treated.data[sig_treated_untreated_up$gene, ])
untreated.data.df <- data.frame(untreated.data[sig_treated_untreated_up$gene, ])

# downsample treated
treated_keep <- sample(ncol(treated.data.df), cells_per_group)
treated.data.DS <- treated.data.df[, treated_keep]
treated.umap.DS <- data.frame(treated.umap[treated_keep, ])
treated.cluster <- data.frame("X" = rep(1, cells_per_group))

# downsample untreated
untreated_keep <- sample(ncol(untreated.data.df), cells_per_group)
untreated.data.DS <- untreated.data.df[, untreated_keep]
untreated.umap.DS <- data.frame(untreated.umap[untreated_keep, ])
untreated.cluster <- data.frame("X" = rep(0, cells_per_group))

# combine
markers.df <- cbind(treated.data.DS, untreated.data.DS)
umap.df <- rbind(treated.umap.DS, untreated.umap.DS)
clusters.df <- rbind(treated.cluster, untreated.cluster)

row.names(clusters.df) <- row.names(umap.df)
colnames(markers.df) <- row.names(umap.df)

# write COMET input files
write.table(markers.df, paste0(figdir3,"treated_vs_untreated_markers.txt"),
            sep="\t", row.names=TRUE, quote=FALSE)

write.table(umap.df, paste0(figdir3,"treated_vs_untreated_umap.txt"),
            sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

write.table(clusters.df, paste0(figdir3,"treated_vs_untreated_clusters.txt"),
            sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

write.table(sig_treated_untreated_up$gene,
            paste0(figdir3,"treated_vs_untreated_genes.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# checks
intersect(colnames(markers.df), rownames(clusters.df))

genes <- intersect(sig_treated_untreated_up$gene, rownames(markers.df))
all <- intersect(rownames(clusters.df), rownames(umap.df))


#########################################################
#### extract data for COMET in early vs advanced ####
#########################################################

figdir4 <- "/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/COMET_EvA/"

Idents(panH) <- "progression_group"

# DEG between advanced and early cells
early_vs_advanced.markers <- FindMarkers(panH, ident.1 = "Advanced", ident.2 = "Early")
early_vs_advanced.markers$gene <- rownames(early_vs_advanced.markers)

# significant genes
sig_early_advanced_genes <- early_vs_advanced.markers %>% 
  filter(p_val_adj <= 0.05)

# upregulated genes in advanced
sig_early_advanced_up <- sig_early_advanced_genes %>% 
  filter(avg_log2FC > 1)

# normalized expression matrix
panH_counts_matrix <- GetAssayData(panH, layer = "data", assay = "SCT")

# advanced cells
advanced_cells <- panH$progression_group == "Advanced"
advanced.data <- panH_counts_matrix[,advanced_cells]

# UMAP
panH_embeddings <- Embeddings(panH, reduction = "umap")
advanced.umap <- panH_embeddings[advanced_cells,]

# early cells
early_cells <- panH$progression_group == "Early"
early.data <- panH_counts_matrix[,early_cells]
early.umap <- panH_embeddings[early_cells,]

# COMET max = 65k cells, balanced between groups
set.seed(1)  # for reproducibility
cells_per_group <- min(32500, ncol(advanced.data), ncol(early.data))

# convert to dataframe
advanced.data.df <- data.frame(advanced.data[sig_early_advanced_up$gene, ])
early.data.df <- data.frame(early.data[sig_early_advanced_up$gene, ])

# downsample advanced
advanced_keep <- sample(ncol(advanced.data.df), cells_per_group)
advanced.data.DS <- advanced.data.df[, advanced_keep]
advanced.umap.DS <- data.frame(advanced.umap[advanced_keep, ])
advanced.cluster <- data.frame("X" = rep(1, cells_per_group))

# downsample early
early_keep <- sample(ncol(early.data.df), cells_per_group)
early.data.DS <- early.data.df[, early_keep]
early.umap.DS <- data.frame(early.umap[early_keep, ])
early.cluster <- data.frame("X" = rep(0, cells_per_group))

# combine
markers.df <- cbind(advanced.data.DS, early.data.DS)
umap.df <- rbind(advanced.umap.DS, early.umap.DS)
clusters.df <- rbind(advanced.cluster, early.cluster)

row.names(clusters.df) <- row.names(umap.df)
colnames(markers.df) <- row.names(umap.df)

# write COMET input files
write.table(markers.df, paste0(figdir4,"early_vs_advanced_markers.txt"),
            sep="\t", row.names=TRUE, quote=FALSE)

write.table(umap.df, paste0(figdir4,"early_vs_advanced_umap.txt"),
            sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

write.table(clusters.df, paste0(figdir4,"early_vs_advanced_clusters.txt"),
            sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

write.table(sig_early_advanced_up$gene,
            paste0(figdir4,"early_vs_advanced_genes.txt"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# checks
intersect(colnames(markers.df), rownames(clusters.df))

genes <- intersect(sig_early_advanced_up$gene, rownames(markers.df))
all <- intersect(rownames(clusters.df), rownames(umap.df))

####################
#####Meta vs Prim######
#######################

# FindMarkers: Metastasis vs Primary
Idents(panH) <- "tissue.type"
metavsprim.markers <- FindMarkers(panH, ident.1 = "Metastasis", ident.2 = "Primary")
metavsprim.markers$gene <- rownames(metavsprim.markers)

# Significant genes
sig_metavsprim_genes <- metavsprim.markers %>% filter(p_val_adj <= 0.05)

# Upregulated genes
sig_metavsprim_path_up <- sig_metavsprim_genes %>% filter(avg_log2FC > 1)

# Convert SYMBOL to ENTREZID
sig_metavsprim_path_up_map <- bitr(
  geneID = sig_metavsprim_path_up$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

# COMET ranked genes
COMET_genes <- read.csv('/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/COMET_MvP/clus2_output/data/cluster_1_singleton_all_ranked.csv')

# Pathways table
msig_Hpath_pvm_upgenes_df <- read.csv("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Prim_v_Meta.csv")

########################
#HALLMARK_ALLOGRAFT_REJECTION
########################
Allo_genes <- msig_Hpath_pvm_upgenes_df$geneID[msig_Hpath_pvm_upgenes_df$Description == "HALLMARK_ALLOGRAFT_REJECTION"]
Allo_genes <- unlist(strsplit(Allo_genes, "/")) %>% trimws()

Allo_geneNames <- sig_metavsprim_path_up_map$ENTREZID %in% Allo_genes
Allo_geneNames <- sig_metavsprim_path_up_map[Allo_geneNames, "SYMBOL"]  

Allo_COMETrank <- COMET_genes$gene_1 %in% Allo_geneNames
Allo_COMETrank <- COMET_genes[Allo_COMETrank, "rank"]

Allo_genes_df <- data.frame("ENTREZID" = Allo_genes, "Gene" = Allo_geneNames)


Allo_genes_final <- Allo_genes_df %>% left_join(
  COMET_genes %>% dplyr::select(Gene = gene_1, COMETrank = rank),
  by = "Gene"
)

########################
#HALLMARK ESTROGEN RESPONSE LATE
########################

ERL_genes <- msig_Hpath_pvm_upgenes_df$geneID[msig_Hpath_pvm_upgenes_df$Description == "HALLMARK_ESTROGEN_RESPONSE_LATE"]
ERL_genes <- unlist(strsplit(ERL_genes, "/")) %>% trimws()

ERL_geneNames <- sig_metavsprim_path_up_map$ENTREZID %in% ERL_genes
ERL_geneNames <- sig_metavsprim_path_up_map[ERL_geneNames, "SYMBOL"]  

ERL_COMETrank <- COMET_genes$gene_1 %in% ERL_geneNames
ERL_COMETrank <- COMET_genes[ERL_COMETrank, "rank"]

ERL_genes_df <- data.frame("ENTREZID" = ERL_genes, "Gene" = ERL_geneNames)


ERL_genes_final <- ERL_genes_df %>% left_join(
  COMET_genes %>% dplyr::select(Gene = gene_1, COMETrank = rank),
  by = "Gene"
)

########################
#HALLMARK E2F TARGETS
########################

E2F_genes <- msig_Hpath_pvm_upgenes_df$geneID[msig_Hpath_pvm_upgenes_df$Description == "HALLMARK_E2F_TARGETS"]
E2F_genes <- unlist(strsplit(E2F_genes, "/")) %>% trimws()

E2F_geneNames <- sig_metavsprim_path_up_map$ENTREZID %in% E2F_genes
E2F_geneNames <- sig_metavsprim_path_up_map[E2F_geneNames, "SYMBOL"]  

E2F_COMETrank <- COMET_genes$gene_1 %in% E2F_geneNames
E2F_COMETrank <- COMET_genes[E2F_COMETrank, "rank"]

E2F_genes_df <- data.frame("ENTREZID" = E2F_genes, "Gene" = E2F_geneNames)


E2F_genes_final <- E2F_genes_df %>% left_join(
  COMET_genes %>% dplyr::select(Gene = gene_1, COMETrank = rank),
  by = "Gene"
)

#Downregulated
# Downregulated genes
sig_metavsprim_path_down <- sig_metavsprim_genes %>% filter(avg_log2FC < -1)

# Convert SYMBOL to ENTREZID
sig_metavsprim_path_down_map <- bitr(
  geneID = sig_metavsprim_path_down$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

# Pathways table
msig_Hpath_mvp_downgenes_df <- read.csv("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Prim_v_Meta.csv")


########################
#HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
########################


EMT_genes <- msig_Hpath_mvp_downgenes_df$geneID[msig_Hpath_mvp_downgenes_df$Description == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
EMT_genes <- unlist(strsplit(EMT_genes, "/")) %>% trimws()

EMT_geneNames <- sig_metavsprim_path_down_map$ENTREZID %in% EMT_genes
EMT_geneNames <- sig_metavsprim_path_down_map[EMT_geneNames, "SYMBOL"]  


EMT_genes_df <- data.frame("ENTREZID" = EMT_genes, "Gene" = EMT_geneNames)


########################
#HALLMARK_ANGIOGENESIS
########################


Ang_genes <- msig_Hpath_mvp_downgenes_df$geneID[msig_Hpath_mvp_downgenes_df$Description == "HALLMARK_ANGIOGENESIS"]
Ang_genes <- unlist(strsplit(Ang_genes, "/")) %>% trimws()

Ang_geneNames <- sig_metavsprim_path_down_map$ENTREZID %in% Ang_genes
Ang_geneNames <- sig_metavsprim_path_down_map[Ang_geneNames, "SYMBOL"]  


Ang_genes_df <- data.frame("ENTREZID" = Ang_genes, "Gene" = Ang_geneNames)


########################
#HALLMARK_MYOGENESIS
########################


Myo_genes <- msig_Hpath_mvp_downgenes_df$geneID[msig_Hpath_mvp_downgenes_df$Description == "HALLMARK_MYOGENESIS"]
Myo_genes <- unlist(strsplit(Myo_genes, "/")) %>% trimws()

Myo_geneNames <- sig_metavsprim_path_down_map$ENTREZID %in% Myo_genes
Myo_geneNames <- sig_metavsprim_path_down_map[Myo_geneNames, "SYMBOL"]  


Myo_genes_df <- data.frame("ENTREZID" = Myo_genes, "Gene" = Myo_geneNames)


####################
#####Treated vs Untreated######
#######################

# FindMarkers: Metastasis vs Primary
Idents(panH) <- "treatment"
tvsu.markers <- FindMarkers(panH, ident.1 = "Treated", ident.2 = "Untreated")
tvsu.markers$gene <- rownames(tvsu.markers)

# Significant genes
sig_tvsu_genes <- tvsu.markers %>% filter(p_val_adj <= 0.05)

# Upregulated genes
sig_tvsu_path_up <- sig_tvsu_genes %>% filter(avg_log2FC > 1)

# Convert SYMBOL to ENTREZID
sig_tvsu_path_up_map <- bitr(
  geneID = sig_tvsu_path_up$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

# COMET ranked genes
COMET_genes <- read.csv('/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/COMET_TvU/treated_output/data/cluster_1_singleton_all_ranked.csv')

# Pathways table
msig_Hpath_tvsu_upgenes_df <- read.csv("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Treated_v_Untreated.csv")

########################
#HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
########################
EMT2_genes <- msig_Hpath_tvsu_upgenes_df$geneID[msig_Hpath_tvsu_upgenes_df$Description == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
EMT2_genes <- unlist(strsplit(EMT2_genes, "/")) %>% trimws()

EMT2_geneNames <- sig_tvsu_path_up_map$ENTREZID %in% EMT2_genes
EMT2_geneNames <- sig_tvsu_path_up_map[EMT2_geneNames, "SYMBOL"]  

EMT2_COMETrank <- COMET_genes$gene_1 %in% EMT2_geneNames
EMT2_COMETrank <- COMET_genes[EMT2_COMETrank, "rank"]

EMT2_genes_df <- data.frame("ENTREZID" = EMT2_genes, "Gene" = EMT2_geneNames)


EMT2_genes_final <- EMT2_genes_df %>% left_join(
  COMET_genes %>% dplyr::select(Gene = gene_1, COMETrank = rank),
  by = "Gene"
)



########################
#HALLMARK_ANGIOGENESIS
########################
Ang2_genes <- msig_Hpath_tvsu_upgenes_df$geneID[msig_Hpath_tvsu_upgenes_df$Description == "HALLMARK_ANGIOGENESIS"]
Ang2_genes <- unlist(strsplit(Ang2_genes, "/")) %>% trimws()

Ang2_geneNames <- sig_tvsu_path_up_map$ENTREZID %in% Ang2_genes
Ang2_geneNames <- sig_tvsu_path_up_map[Ang2_geneNames, "SYMBOL"]  

Ang2_COMETrank <- COMET_genes$gene_1 %in% Ang2_geneNames
Ang2_COMETrank <- COMET_genes[Ang2_COMETrank, "rank"]

Ang2_genes_df <- data.frame("ENTREZID" = Ang2_genes, "Gene" = Ang2_geneNames)

Ang2_genes_final <- Ang2_genes_df %>% left_join(
  COMET_genes %>% dplyr::select(Gene = gene_1, COMETrank = rank),
  by = "Gene"
)

########################
#HALLMARK_MYOGENESIS
########################
Myo2_genes <- msig_Hpath_tvsu_upgenes_df$geneID[msig_Hpath_tvsu_upgenes_df$Description == "HALLMARK_MYOGENESIS"]
Myo2_genes <- unlist(strsplit(Myo2_genes, "/")) %>% trimws()

Myo2_geneNames <- sig_tvsu_path_up_map$ENTREZID %in% Myo2_genes
Myo2_geneNames <- sig_tvsu_path_up_map[Myo2_geneNames, "SYMBOL"]  

Myo2_COMETrank <- COMET_genes$gene_1 %in% Myo2_geneNames
Myo2_COMETrank <- COMET_genes[Myo2_COMETrank, "rank"]

Myo2_genes_df <- data.frame("ENTREZID" = Myo2_genes, "Gene" = Myo2_geneNames)

Myo2_genes_final <- Myo2_genes_df %>% left_join(
  COMET_genes %>% dplyr::select(Gene = gene_1, COMETrank = rank),
  by = "Gene"
)

#Downregulated
# Downregulated genes
sig_tvsu_path_down <- sig_tvsu_genes %>% filter(avg_log2FC < -1)

# Convert SYMBOL to ENTREZID
sig_tvsu_path_down_map <- bitr(
  geneID = sig_tvsu_path_down$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

# Pathways table
msig_Hpath_tvu_downgenes_df <- read.csv("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Treated_v_Untreated.csv")


########################
#HALLMARK_ALLOGRAFT_REJECTION
########################


Allo2_genes <- msig_Hpath_tvu_downgenes_df$geneID[msig_Hpath_tvu_downgenes_df$Description == "HALLMARK_ALLOGRAFT_REJECTION"]
Allo2_genes <- unlist(strsplit(Allo2_genes, "/")) %>% trimws()

Allo2_geneNames <- sig_tvsu_path_down_map$ENTREZID %in% Allo2_genes
Allo2_geneNames <- sig_tvsu_path_down_map[Allo2_geneNames, "SYMBOL"]  


Allo2_genes_df <- data.frame("ENTREZID" = Allo2_genes, "Gene" = Allo2_geneNames)

########################
#HALLMARK_COMPLEMENT
########################


Comp_genes <- msig_Hpath_tvu_downgenes_df$geneID[msig_Hpath_tvu_downgenes_df$Description == "HALLMARK_COMPLEMENT"]
Comp_genes <- unlist(strsplit(Comp_genes, "/")) %>% trimws()

Comp_geneNames <- sig_tvsu_path_down_map$ENTREZID %in% Comp_genes
Comp_geneNames <- sig_tvsu_path_down_map[Comp_geneNames, "SYMBOL"]  


Comp_genes_df <- data.frame("ENTREZID" = Comp_genes, "Gene" = Comp_geneNames)


########################
#HALLMARK_INFLAMMATORY_RESPONSE
########################


Inf_genes <- msig_Hpath_tvu_downgenes_df$geneID[msig_Hpath_tvu_downgenes_df$Description == "HALLMARK_INFLAMMATORY_RESPONSE"]
Inf_genes <- unlist(strsplit(Inf_genes, "/")) %>% trimws()

Inf_geneNames <- sig_tvsu_path_down_map$ENTREZID %in% Inf_genes
Inf_geneNames <- sig_tvsu_path_down_map[Inf_geneNames, "SYMBOL"]  


Inf_genes_df <- data.frame("ENTREZID" = Inf_genes, "Gene" = Inf_geneNames)

####################
####Advanced vs Early######
#######################

# FindMarkers: Advanced vs Early
Idents(panH) <- "progression_group"
avse.markers <- FindMarkers(panH, ident.1 = "Advanced", ident.2 = "Early")
avse.markers$gene <- rownames(avse.markers)

# Significant genes
sig_avse_genes <- avse.markers %>% filter(p_val_adj <= 0.05)

# Upregulated genes
sig_avse_path_up <- sig_avse_genes %>% filter(avg_log2FC > 1)

# Convert SYMBOL to ENTREZID
sig_avse_path_up_map <- bitr(
  geneID = sig_avse_path_up$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

# COMET ranked genes
COMET_genes <- read.csv('/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/COMET_EvA/early_vs_advanced_output/data/cluster_1_singleton_all_ranked.csv')

# Pathways table
msig_Hpath_avse_upgenes_df <- read.csv("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/mSigHpathways-DEG0.05PvalueANDlower1FC-Downregulated_Higher_Early_v_Advanced.csv")

########################
#HALLMARK_GLYCOLYSIS
########################
Gly_genes <- msig_Hpath_avse_upgenes_df$geneID[msig_Hpath_avse_upgenes_df$Description == "HALLMARK_GLYCOLYSIS"]
Gly_genes <- unlist(strsplit(Gly_genes, "/")) %>% trimws()

Gly_geneNames <- sig_avse_path_up_map$ENTREZID %in% Gly_genes
Gly_geneNames <- sig_avse_path_up_map[Gly_geneNames, "SYMBOL"]  

Gly_COMETrank <- COMET_genes$gene_1 %in% Gly_geneNames
Gly_COMETrank <- COMET_genes[Gly_COMETrank, "rank"]

Gly_genes_df <- data.frame("ENTREZID" = Gly_genes, "Gene" = Gly_geneNames)


Gly_genes_final <- Gly_genes_df %>% left_join(
  COMET_genes %>% dplyr::select(Gene = gene_1, COMETrank = rank),
  by = "Gene"
)

########################
#HALLMARK_INFLAMMATORY_RESPONSE
########################
Inf3_genes <- msig_Hpath_avse_upgenes_df$geneID[msig_Hpath_avse_upgenes_df$Description == "HALLMARK_INFLAMMATORY_RESPONSE"]
Inf3_genes <- unlist(strsplit(Inf3_genes, "/")) %>% trimws()

Inf3_geneNames <- sig_avse_path_up_map$ENTREZID %in% Inf3_genes
Inf3_geneNames <- sig_avse_path_up_map[Inf3_geneNames, "SYMBOL"]  

Inf3_COMETrank <- COMET_genes$gene_1 %in% Inf3_geneNames
Inf3_COMETrank <- COMET_genes[Inf3_COMETrank, "rank"]

Inf3_genes_df <- data.frame("ENTREZID" = Inf3_genes, "Gene" = Inf3_geneNames)


Inf3_genes_final <- Inf3_genes_df %>% left_join(
  COMET_genes %>% dplyr::select(Gene = gene_1, COMETrank = rank),
  by = "Gene"
)

########################
#HALLMARK_ESTROGEN_RESPONSE_LATE
########################
ERL3_genes <- msig_Hpath_avse_upgenes_df$geneID[msig_Hpath_avse_upgenes_df$Description == "HALLMARK_ESTROGEN_RESPONSE_LATE"]
ERL3_genes <- unlist(strsplit(ERL3_genes, "/")) %>% trimws()

ERL3_geneNames <- sig_avse_path_up_map$ENTREZID %in% ERL3_genes
ERL3_geneNames <- sig_avse_path_up_map[ERL3_geneNames, "SYMBOL"]  

ERL3_COMETrank <- COMET_genes$gene_1 %in% ERL3_geneNames
ERL3_COMETrank <- COMET_genes[ERL3_COMETrank, "rank"]

ERL3_genes_df <- data.frame("ENTREZID" = ERL3_genes, "Gene" = ERL3_geneNames)


ERL3_genes_final <- ERL3_genes_df %>% left_join(
  COMET_genes %>% dplyr::select(Gene = gene_1, COMETrank = rank),
  by = "Gene"
)

#Downregulated
# Downregulated genes
sig_avse_path_down <- sig_avse_genes %>% filter(avg_log2FC < -1)

# Convert SYMBOL to ENTREZID
sig_avse_path_down_map <- bitr(
  geneID = sig_avse_path_down$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

# Pathways table
msig_Hpath_ave_downgenes_df <- read.csv("/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/mSigHpathways-DEG0.05PvalueANDhigher1FC-Upregulated_Higher_Early_v_Advanced.csv")


########################
#HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
########################


EMT3_genes <- msig_Hpath_ave_downgenes_df$geneID[msig_Hpath_ave_downgenes_df$Description == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
EMT3_genes <- unlist(strsplit(EMT3_genes, "/")) %>% trimws()

EMT3_geneNames <- sig_avse_path_down_map$ENTREZID %in% EMT3_genes
EMT3_geneNames <- sig_avse_path_down_map[EMT3_geneNames, "SYMBOL"]  


EMT3_genes_df <- data.frame("ENTREZID" = EMT3_genes, "Gene" = EMT3_geneNames)


########################
#HALLMARK_ANGIOGENESIS
########################


Ang3_genes <- msig_Hpath_ave_downgenes_df$geneID[msig_Hpath_ave_downgenes_df$Description == "HALLMARK_ANGIOGENESIS"]
Ang3_genes <- unlist(strsplit(Ang3_genes, "/")) %>% trimws()

Ang3_geneNames <- sig_avse_path_down_map$ENTREZID %in% Ang3_genes
Ang3_geneNames <- sig_avse_path_down_map[Ang3_geneNames, "SYMBOL"]  


Ang3_genes_df <- data.frame("ENTREZID" = Ang3_genes, "Gene" = Ang3_geneNames)


########################
#HALLMARK_MYOGENESIS
########################


Myo3_genes <- msig_Hpath_ave_downgenes_df$geneID[msig_Hpath_ave_downgenes_df$Description == "HALLMARK_MYOGENESIS"]
Myo3_genes <- unlist(strsplit(Myo3_genes, "/")) %>% trimws()

Myo3_geneNames <- sig_avse_path_down_map$ENTREZID %in% Myo3_genes
Myo3_geneNames <- sig_avse_path_down_map[Myo3_geneNames, "SYMBOL"]  


Myo3_genes_df <- data.frame("ENTREZID" = Myo3_genes, "Gene" = Myo3_geneNames)


#PLots for COMET genes
#MvP

# Create the plot and store in a variable
p <- DimPlot(panH, reduction = "umap", group.by = "tissue.type")

# Save to file
ggsave(filename = paste0(figdir, "met_vs_primary_umap.png"), 
       plot = p, width = 6, height = 5, dpi = 300)

p <- FeaturePlot(panH, features = "RAMP1", reduction = "umap")
ggsave(filename = paste0(figdir, "met_vs_primary_RAMP1_feature.png"), 
       plot = p, width = 6, height = 5, dpi = 300)

#AvE

# Create the plot and store in a variable
p <- DimPlot(panH, reduction = "umap", group.by = "progression_group")

# Save to file
ggsave(filename = paste0(figdir, "adv_vs_early_umap.png"), 
       plot = p, width = 6, height = 5, dpi = 300)

p <- FeaturePlot(panH, features = "DDX3Y", reduction = "umap")
ggsave(filename = paste0(figdir, "adv_vs_early__DDX3Y_feature.png"), 
       plot = p, width = 6, height = 5, dpi = 300)

p <- FeaturePlot(panH, features = "RPS4Y1", reduction = "umap")
ggsave(filename = paste0(figdir, "adv_vs_early_RPS4Y1_feature.png"), 
       plot = p, width = 6, height = 5, dpi = 300)

#TvU

# Create the plot and store in a variable
p <- DimPlot(panH, reduction = "umap", group.by = "treatment")

# Save to file
ggsave(filename = paste0(figdir, "treated_vs_untreated_umap.png"), 
       plot = p, width = 6, height = 5, dpi = 300)

p <- FeaturePlot(panH, features = "MGP", reduction = "umap")
ggsave(filename = paste0(figdir, "treated_vs_untreated_umap_MGP.png"), 
       plot = p, width = 6, height = 5, dpi = 300)



#Make Chord plots:
# Make Chord plots

library(circlize)
library(dplyr)

figdir <- "/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Data
# -----------------------------
txt <- "Gene\tAllo\tERL\tE2F\tEMT\tAng\tMyo\tlogFC
PRF1\t1\t0\t0\t0\t0\t0\t1.349565
CSK\t1\t0\t0\t0\t0\t0\t1.039931
PTPN6\t1\t1\t0\t0\t0\t0\t1.115492
PLAAT3\t0\t1\t0\t0\t0\t0\t1.309586
MDK\t0\t1\t0\t0\t0\t0\t1.172027
HMGA1\t0\t0\t1\t0\t0\t0\t1.40499
NME1\t1\t0\t1\t0\t0\t0\t1.040377
TK1\t0\t0\t1\t0\t0\t0\t1.438437
COL1A2\t0\t0\t0\t1\t0\t0\t-3.929201
MGP\t0\t0\t0\t1\t0\t0\t-3.749536
COL3A1\t0\t0\t0\t1\t1\t1\t-3.923201
LUM\t0\t0\t0\t1\t1\t0\t-4.468499
VCAN\t0\t0\t0\t1\t1\t0\t-1.054982
SPARC\t0\t0\t0\t1\t0\t1\t-3.013571
COL1A1\t0\t0\t0\t1\t0\t1\t-3.59101"

df <- read.delim(text = txt, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

pathways <- c("Allo", "ERL", "E2F", "EMT", "Ang", "Myo")

# -----------------------------
# Build edges
# -----------------------------
edges <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
  gene <- df$Gene[i]
  fc <- df$logFC[i]
  present <- pathways[as.numeric(df[i, pathways]) == 1]
  
  if (length(present) == 0) return(NULL)
  
  data.frame(
    from = gene,
    to = present,
    value = abs(fc),
    logFC = fc,
    stringsAsFactors = FALSE
  )
}))

# Order: genes first, then pathways
gene_order <- df %>% arrange(logFC) %>% pull(Gene)
sector_order <- c(gene_order, pathways)

# -----------------------------
# Colors
# -----------------------------
grid.col <- c(
  Allo = "#8C2D04",   
  ERL  = "#E6550D",   
  E2F  = "#FDBE85",
  EMT  = "#B7E4C7",   
  Ang  = "#52B788",   
  Myo  = "#1B4332"
)

# Continuous gene colors based on logFC
gene_cols <- setNames(
  col_fun(df$logFC),
  df$Gene
)

grid.col <- c(gene_cols, grid.col)

# Link colors follow pathway color
link_cols <- grid.col[edges$to]

# -----------------------------
# Continuous logFC color function
# -----------------------------
lim <- ceiling(max(abs(df$logFC)))
col_fun <- colorRamp2(
  c(-lim, 0, lim),
  c("#3A7DCE", "white", "#D64550")
)

# -----------------------------
# Save figure
# -----------------------------
pdf(file.path(figdir, "chord_with_logFC.pdf"), width = 10, height = 12)

layout(matrix(c(1, 2), nrow = 2), heights = c(4.5, 1))
par(oma = c(0, 0, 0, 0))

# -----------------------------
# CHORD PLOT
# -----------------------------
par(mar = c(1, 1, 2, 1))

circos.clear()

chordDiagram(
  x = edges[, c("from", "to", "value")],
  order = sector_order,
  grid.col = grid.col,
  col = link_cols,
  transparency = 0.15,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.08),
  directional = 0
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    circos.text(
      x = mean(xlim),
      y = ylim[1] + 0.1,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.75
    )
  },
  bg.border = NA
)

title("Gene–Program Chord Plot", cex.main = 1.2)

# -----------------------------
# LOGFC GRADIENT LEGEND
# -----------------------------
par(mar = c(2, 3, 2, 3))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

# gradient strip placement
x_left  <- 0.22
x_right <- 0.85
y_bot   <- 0.42
y_top   <- 0.72

vals <- seq(-lim, lim, length.out = 300)
cols <- col_fun(vals)

# draw smooth gradient
for (i in seq_len(length(vals) - 1)) {
  x1 <- x_left + (i - 1) / (length(vals) - 1) * (x_right - x_left)
  x2 <- x_left + i / (length(vals) - 1) * (x_right - x_left)
  rect(x1, y_bot, x2, y_top, col = cols[i], border = NA)
}

# legend label
text(0.08, (y_bot + y_top) / 2, "logFC", cex = 1.4, xpd = NA)

# ticks
tick_vals <- pretty(c(-lim, lim), n = 5)
tick_pos <- x_left + (tick_vals + lim) / (2 * lim) * (x_right - x_left)

segments(tick_pos, y_bot - 0.03, tick_pos, y_bot, lwd = 1)
text(tick_pos, y_bot - 0.08, labels = tick_vals, cex = 1)

dev.off()
circos.clear()


# Make Chord plots

library(circlize)
library(dplyr)

figdir <- "/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Data
# -----------------------------
txt <- "Gene\tEMT\tAng\tMyo\tAllo\tComp\tInf2\tlogFC
MGP\t1\t0\t0\t0\t0\t0\t2.445379
FBLN1\t1\t0\t0\t0\t0\t0\t3.101355
DCN\t1\t0\t0\t0\t0\t0\t2.623718
LUM\t1\t1\t0\t0\t0\t0\t1.92608
COL3A1\t1\t1\t1\t0\t0\t0\t2.216002
COL5A2\t1\t1\t0\t0\t0\t0\t2.005429
COL1A1\t1\t0\t1\t0\t0\t0\t1.799471
SPARC\t1\t0\t1\t0\t0\t0\t1.676179
GZMA\t0\t0\t0\t1\t1\t0\t-1.364849
ITGB2\t0\t0\t0\t1\t0\t0\t-1.082428
CCL4\t0\t0\t0\t1\t0\t0\t-1.168113
GZMK\t0\t0\t0\t0\t1\t0\t-1.461696
APOC1\t0\t0\t0\t0\t1\t0\t-1.504053
IFITM1\t0\t0\t0\t0\t0\t1\t-1.11259
LCK\t0\t0\t0\t1\t1\t1\t-1.049508
IL2RB\t0\t0\t0\t1\t0\t1\t-1.30021"

df <- read.delim(text = txt, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

pathways <- c("EMT", "Ang", "Myo", "Allo", "Comp", "Inf2")

# -----------------------------
# Build edges
# -----------------------------
edges <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
  gene <- df$Gene[i]
  fc <- df$logFC[i]
  present <- pathways[as.numeric(df[i, pathways]) == 1]
  
  if (length(present) == 0) return(NULL)
  
  data.frame(
    from = gene,
    to = present,
    value = abs(fc),
    logFC = fc,
    stringsAsFactors = FALSE
  )
}))

# Order
gene_order <- df %>% arrange(logFC) %>% pull(Gene)
sector_order <- c(gene_order, pathways)

# -----------------------------
# Colors
# -----------------------------

# Pathway colors
pathway_cols <- c(
  EMT = "#8C2D04",   # dark burnt orange
  Ang = "#E6550D",   # strong orange
  Myo = "#FDBE85",    # light peach
  Allo = "#B7E4C7",   # light green
  Comp = "#52B788",   # medium green
  Inf2 = "#1B4332"   # dark green
)

lim <- max(abs(df$logFC))
col_fun <- colorRamp2(
  c(-lim, 0, lim),
  c("#3A7DCE", "white", "#D64550")
)

# Gene colors follow actual logFC
gene_cols <- setNames(
  col_fun(df$logFC),
  df$Gene
)

# Combined sector colors
grid.col <- c(gene_cols, pathway_cols)

# Link colors follow pathway colors
link_cols <- adjustcolor(pathway_cols[edges$to], alpha.f = 0.7)

# -----------------------------
# Save figure
# -----------------------------
pdf(file.path(figdir, "chord_with_logFC_set2.pdf"), width = 10, height = 12)

layout(matrix(c(1, 2), nrow = 2), heights = c(4.5, 1))
par(oma = c(0, 0, 0, 0))

# -----------------------------
# CHORD PLOT
# -----------------------------
par(mar = c(1, 1, 2, 1))

circos.clear()

chordDiagram(
  x = edges[, c("from", "to", "value")],
  order = sector_order,
  grid.col = grid.col,
  col = link_cols,
  transparency = 0.15,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.08),
  directional = 0
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    circos.text(
      x = mean(xlim),
      y = ylim[1] + 0.1,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.75
    )
  },
  bg.border = NA
)

title("Gene–Program Chord Plot", cex.main = 1.2)

# -----------------------------
# LOGFC GRADIENT LEGEND
# -----------------------------
par(mar = c(2, 3, 2, 3))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

x_left  <- 0.22
x_right <- 0.85
y_bot   <- 0.42
y_top   <- 0.72

vals <- seq(-lim, lim, length.out = 300)
cols <- col_fun(vals)

for (i in seq_len(length(vals) - 1)) {
  x1 <- x_left + (i - 1) / (length(vals) - 1) * (x_right - x_left)
  x2 <- x_left + i / (length(vals) - 1) * (x_right - x_left)
  rect(x1, y_bot, x2, y_top, col = cols[i], border = NA)
}

text(0.08, (y_bot + y_top) / 2, "logFC", cex = 1.4, xpd = NA)

tick_vals <- pretty(c(-lim, lim), n = 5)
tick_pos <- x_left + (tick_vals + lim) / (2 * lim) * (x_right - x_left)

segments(tick_pos, y_bot - 0.03, tick_pos, y_bot, lwd = 1)
text(tick_pos, y_bot - 0.08, labels = tick_vals, cex = 1)

dev.off()
circos.clear()


# Make Chord plots

library(circlize)
library(dplyr)

figdir <- "/blue/ferrallm/antonihalka/Pancreas/scripts/Spring2026/plots/"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Data
# -----------------------------
txt <- "Gene\tGly\tInf3\tERL\tEMT\tAng\tMyo\tlogFC
MET\t1\t1\t0\t0\t0\t0\t1.603652
ERO1A\t1\t0\t0\t0\t0\t0\t1.329461
P4HA1\t1\t0\t0\t0\t0\t0\t1.040996
IFITM1\t0\t1\t0\t0\t0\t0\t1.02775
DCBLD2\t0\t1\t0\t0\t0\t0\t2.104674
PERP\t0\t0\t1\t0\t0\t0\t1.017527
KRT19\t0\t0\t1\t0\t0\t0\t1.148093
LAMC2\t0\t0\t1\t0\t0\t0\t1.141709
COL1A2\t0\t0\t0\t1\t0\t0\t-2.352466
MGP\t0\t0\t0\t1\t0\t0\t-2.018055
COL3A1\t0\t0\t0\t1\t1\t1\t-2.18145
LUM\t0\t0\t0\t1\t1\t0\t-2.147903
FSTL1\t0\t0\t0\t1\t1\t0\t-1.788678
COL1A1\t0\t0\t0\t1\t0\t1\t-2.289156
SPARC\t0\t0\t0\t1\t0\t1\t-1.909015"

df <- read.delim(text = txt, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- trimws(colnames(df))

pathways <- c("Gly", "Inf3", "ERL", "EMT", "Ang", "Myo")

# -----------------------------
# Build edges
# -----------------------------
edges <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
  gene <- df$Gene[i]
  fc <- df$logFC[i]
  present <- pathways[as.numeric(unlist(df[i, pathways])) == 1]
  
  if (length(present) == 0) return(NULL)
  
  data.frame(
    from = gene,
    to = present,
    value = abs(fc),
    logFC = fc,
    stringsAsFactors = FALSE
  )
}))

# Order
gene_order <- df %>% arrange(logFC) %>% pull(Gene)
sector_order <- c(gene_order, pathways)

# -----------------------------
# Colors
# -----------------------------
pathway_cols <- c(
  Gly = "#8C2D04",   # dark burnt orange
  Inf3 = "#E6550D",   # strong orange
  ERL = "#FDBE85",   # light peach
  EMT = "#B7E4C7",   # light green
  Ang = "#52B788",   # medium green
  Myo = "#1B4332"    # dark green
)

lim <- max(abs(df$logFC))
col_fun <- colorRamp2(
  c(-lim, 0, lim),
  c("#3A7DCE", "white", "#D64550")
)

gene_cols <- setNames(
  col_fun(df$logFC),
  df$Gene
)

grid.col <- c(gene_cols, pathway_cols)
link_cols <- adjustcolor(pathway_cols[edges$to], alpha.f = 0.7)

# -----------------------------
# Save figure
# -----------------------------
pdf(file.path(figdir, "chord_with_logFC_set3.pdf"), width = 10, height = 12)

layout(matrix(c(1, 2), nrow = 2), heights = c(4.5, 1))
par(oma = c(0, 0, 0, 0))

# -----------------------------
# CHORD PLOT
# -----------------------------
par(mar = c(1, 1, 2, 1))

circos.clear()

chordDiagram(
  x = edges[, c("from", "to", "value")],
  order = sector_order,
  grid.col = grid.col,
  col = link_cols,
  transparency = 0.15,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.08),
  directional = 0
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    circos.text(
      x = mean(xlim),
      y = ylim[1] + 0.1,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.75
    )
  },
  bg.border = NA
)

title("Gene–Program Chord Plot", cex.main = 1.2)

# -----------------------------
# LOGFC GRADIENT LEGEND
# -----------------------------
par(mar = c(2, 3, 2, 3))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

x_left  <- 0.22
x_right <- 0.85
y_bot   <- 0.42
y_top   <- 0.72

vals <- seq(-lim, lim, length.out = 300)
cols <- col_fun(vals)

for (i in seq_len(length(vals) - 1)) {
  x1 <- x_left + (i - 1) / (length(vals) - 1) * (x_right - x_left)
  x2 <- x_left + i / (length(vals) - 1) * (x_right - x_left)
  rect(x1, y_bot, x2, y_top, col = cols[i], border = NA)
}

text(0.08, (y_bot + y_top) / 2, "logFC", cex = 1.4, xpd = NA)

tick_vals <- pretty(c(-lim, lim), n = 5)
tick_pos <- x_left + (tick_vals + lim) / (2 * lim) * (x_right - x_left)

segments(tick_pos, y_bot - 0.03, tick_pos, y_bot, lwd = 1)
text(tick_pos, y_bot - 0.08, labels = tick_vals, cex = 1)

dev.off()
circos.clear()

