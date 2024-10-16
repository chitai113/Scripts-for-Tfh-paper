library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
library(HGNChelper)
library(limma)
library(SeuratData)
library(openxlsx)
library(tidyverse)
library(writexl)

# Import the data 
A4_RA.data <- Read10X(data.dir = "/media/user/Elements/OSU_A3_TO_A6/run_count_A4_RA/outs/filtered_feature_bc_matrix")
A3_HC.data <- Read10X(data.dir = "/media/user/Elements/OSU_A3_TO_A6/run_count_A3_HC/outs/filtered_feature_bc_matrix")
A6_RA.data <- Read10X(data.dir = "/media/user/Elements/OSU_A3_TO_A6/run_count_A6_RA/outs/filtered_feature_bc_matrix")
A5_HC.data <- Read10X(data.dir = "/media/user/Elements/OSU_A3_TO_A6/run_count_A5_HC/outs/filtered_feature_bc_matrix")

A2_RA.data <- Read10X(data.dir = "/media/user/Elements/ravshc/RA/outs")
A1_HC.data <- Read10X(data.dir = "/media/user/Elements/ravshc/HC/outs")


#======================================Processing RA data=============================================================
#=====================================================================================================================

# Processing A2 RA
A2_RA<-CreateSeuratObject(counts = A2_RA.data,project = "RA",min.cells = 3,min.features = 200)
A2_RA$Condition="RA"
A2_RA[["percent.mt"]]<-PercentageFeatureSet(A2_RA,pattern = "mt-")
A2_RA<-subset(A2_RA,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
A2_RA<-NormalizeData(A2_RA, normalization.method = "LogNormalize", scale.factor = 10000)
A2_RA<-FindVariableFeatures(A2_RA,selection.method = "vst",nfeatures = 2000)

all.genes <- rownames(A2_RA)
A2_RA <- ScaleData(A2_RA, features = all.genes)

A2_RA <- RunPCA(A2_RA, features = VariableFeatures(object = A2_RA))
A2_RA <- FindNeighbors(A2_RA, dims = 1:10)
A2_RA <- FindClusters(A2_RA, resolution = 0.5)
A2_RA <- RunUMAP(A2_RA, dims = 1:10)

saveRDS(A2_RA,file="A2_RA") #Save processed A2 RA data

#Assign Idents
Idents(A2_RA) <- "RA"
head(A2_RA)





# Processing A4 RA
A4_RA<-CreateSeuratObject(counts = A4_RA.data,project = "RA",min.cells = 3,min.features = 200)
A4_RA$Condition="RA"
A4_RA[["percent.mt"]]<-PercentageFeatureSet(A4_RA,pattern = "mt-")
A4_RA<-subset(A4_RA,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
A4_RA<-NormalizeData(A4_RA, normalization.method = "LogNormalize", scale.factor = 10000)
A4_RA<-FindVariableFeatures(A4_RA,selection.method = "vst",nfeatures = 2000)

all.genes <- rownames(A4_RA)
A4_RA <- ScaleData(A4_RA, features = all.genes)

A4_RA <- RunPCA(A4_RA, features = VariableFeatures(object = A4_RA))
A4_RA <- FindNeighbors(A4_RA, dims = 1:10)
A4_RA <- FindClusters(A4_RA, resolution = 0.5)
A4_RA <- RunUMAP(A4_RA, dims = 1:10)

saveRDS(A4_RA,file="A4_RA")

#Assign Idents
Idents(A4_RA) <- "RA"
head(A4_RA)




# Processing A6 RA
A6_RA<-CreateSeuratObject(counts = A6_RA.data,project = "RA",min.cells = 3,min.features = 200)
A6_RA$Condition="RA"
A6_RA[["percent.mt"]]<-PercentageFeatureSet(A6_RA,pattern = "mt-")
A6_RA<-subset(A6_RA,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
A6_RA<-NormalizeData(A6_RA, normalization.method = "LogNormalize", scale.factor = 10000)
A6_RA<-FindVariableFeatures(A6_RA,selection.method = "vst",nfeatures = 2000)

all.genes <- rownames(A6_RA)
A6_RA <- ScaleData(A6_RA, features = all.genes)

A6_RA <- RunPCA(A6_RA, features = VariableFeatures(object = A6_RA))
A6_RA <- FindNeighbors(A6_RA, dims = 1:10)
A6_RA <- FindClusters(A6_RA, resolution = 0.5)
A6_RA <- RunUMAP(A6_RA, dims = 1:10)

saveRDS(A6_RA,file="A6_RA")

#Assign Idents
Idents(A6_RA) <- "RA"
head(A6_RA)




#============================= Merge RA A2 to A6 =============================

RA2 <- readRDS("A2_RA")
RA4 <- readRDS("A4_RA")
RA6 <- readRDS("A6_RA")
RA2$Condition="RA2"
RA4$Condition="RA4"
RA6$Condition="RA6"

RA2$Type = "RA"
RA4$Type = "RA"
RA6$Type = "RA"
obj.list<-list(RA2,RA4,RA6)
data.anchors<-FindIntegrationAnchors(obj.list, dims = 1:20)
M.dat<-IntegrateData(anchorset = data.anchors, dims = 1:20)
M.dat <- JoinLayers(M.dat,assay="RNA",layer="counts")
counts <- GetAssayData(M.dat, assay = "RNA",layer = "counts")

# Process merged raw data
RAdat <- CreateSeuratObject(counts = counts,meta.data = M.dat@meta.data)
RAdat[["percent.mt"]]<- PercentageFeatureSet(RAdat, pattern = "^MT-")
RAdat <- subset(RAdat, subset= nFeature_RNA >200 & nFeature_RNA < 4000 & percent.mt <10)
RAdat <- NormalizeData(RAdat)
RAdat <- FindVariableFeatures(RAdat, selection.method = "vst", nFeature = 2000)
all.genes <- rownames(RAdat)
RAdat <- ScaleData(RAdat, verbose = FALSE)
RAdat <- RunPCA(RAdat, npcs = 20, verbose = FALSE)
RAdat<- RunUMAP(RAdat, reduction = "pca", dims = 1:20)
RAdat <- FindNeighbors(RAdat, reduction = "pca", dims = 1:20)
RAdat <- FindClusters(RAdat, resolution = 0.4)

table(RAdat$Condition)
Idents(RAdat) <- "Condition"
saveRDS(RAdat,file = "RA_Merged")









#======================================Processing HC data=============================================================
#=====================================================================================================================

# Processing A1 HC
A1_HC<-CreateSeuratObject(counts = A1_HC.data,project = "HC",min.cells = 3,min.features = 200)
A1_HC$Condition="HC"
A1_HC[["percent.mt"]]<-PercentageFeatureSet(A1_HC,pattern = "mt-")
A1_HC<-subset(A1_HC,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
A1_HC<-NormalizeData(A1_HC, normalization.method = "LogNormalize", scale.factor = 10000)
A1_HC<-FindVariableFeatures(A1_HC,selection.method = "vst",nfeatures = 2000)

all.genes <- rownames(A1_HC)
A1_HC <- ScaleData(A1_HC, features = all.genes)

A1_HC <- RunPCA(A1_HC, features = VariableFeatures(object = A1_HC))
A1_HC <- FindNeighbors(A1_HC, dims = 1:10)
A1_HC <- FindClusters(A1_HC, resolution = 0.5)
A1_HC <- RunUMAP(A1_HC, dims = 1:10)

saveRDS(A1_HC,file="A1_HC") #Save processed A2 RA data

#Assign Idents
Idents(A1_HC) <- "HC"
head(A1_HC)




# Processing A3 HC
A3_HC<-CreateSeuratObject(counts = A3_HC.data,project = "HC",min.cells = 3,min.features = 200)
A3_HC$Condition="HC"
A3_HC[["percent.mt"]]<-PercentageFeatureSet(A3_HC,pattern = "mt-")
A3_HC<-subset(A3_HC,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
A3_HC<-NormalizeData(A3_HC, normalization.method = "LogNormalize", scale.factor = 10000)
A3_HC<-FindVariableFeatures(A3_HC,selection.method = "vst",nfeatures = 2000)

all.genes <- rownames(A3_HC)
A3_HC <- ScaleData(A3_HC, features = all.genes)

A3_HC <- RunPCA(A3_HC, features = VariableFeatures(object = A3_HC))
A3_HC <- FindNeighbors(A3_HC, dims = 1:10)
A3_HC <- FindClusters(A3_HC, resolution = 0.5)
A3_HC <- RunUMAP(A3_HC, dims = 1:10)

saveRDS(A3_HC,file="A3_HC") #Save processed A2 RA data

#Assign Idents
Idents(A3_HC) <- "HC"
head(A3_HC)




# Processing A5 HC
A5_HC<-CreateSeuratObject(counts = A5_HC.data,project = "HC",min.cells = 3,min.features = 200)
A5_HC$Condition="HC"
A5_HC[["percent.mt"]]<-PercentageFeatureSet(A5_HC,pattern = "mt-")
A5_HC<-subset(A5_HC,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
A5_HC<-NormalizeData(A5_HC, normalization.method = "LogNormalize", scale.factor = 10000)
A5_HC<-FindVariableFeatures(A5_HC,selection.method = "vst",nfeatures = 2000)

all.genes <- rownames(A5_HC)
A5_HC <- ScaleData(A5_HC, features = all.genes)

A5_HC <- RunPCA(A5_HC, features = VariableFeatures(object = A5_HC))
A5_HC <- FindNeighbors(A5_HC, dims = 1:10)
A5_HC <- FindClusters(A5_HC, resolution = 0.5)
A5_HC <- RunUMAP(A5_HC, dims = 1:10)

saveRDS(A5_HC,file="A5_HC") #Save processed A2 RA data

#Assign Idents
Idents(A5_HC) <- "HC"
head(A5_HC)




#============================= Merge HC A1 to A5 ============================= 

HC1 <- readRDS("A1_HC")
HC3 <- readRDS("A3_HC")
HC5 <- readRDS("A5_HC")
HC1$Condition="HC1"
HC3$Condition="HC3"
HC5$Condition="HC5"

HC1$Type = "HC"
HC3$Type = "HC"
HC5$Type = "HC"
obj.list<-list(HC1,HC3,HC5)
data.anchors<-FindIntegrationAnchors(obj.list, dims = 1:20)
M.dat<-IntegrateData(anchorset = data.anchors, dims = 1:20)
M.dat <- JoinLayers(M.dat,assay="RNA",layer="counts")
counts <- GetAssayData(M.dat, assay = "RNA",layer = "counts")

HCdat <- CreateSeuratObject(counts = counts,meta.data = M.dat@meta.data)
table(HCdat$Condition)
table(HCdat$Type)
HCdat[["percent.mt"]]<- PercentageFeatureSet(HCdat, pattern = "^MT-")
HCdat <- subset(HCdat, subset= nFeature_RNA >200 & nFeature_RNA < 4000 & percent.mt <10)
HCdat <- NormalizeData(HCdat)
HCdat <- FindVariableFeatures(HCdat, selection.method = "vst", nFeature = 2000)
all.genes <- rownames(HCdat)
HCdat <- ScaleData(HCdat, verbose = FALSE)
HCdat <- RunPCA(HCdat, npcs = 20, verbose = FALSE)
HCdat<- RunUMAP(HCdat, reduction = "pca", dims = 1:20)
HCdat <- FindNeighbors(HCdat, reduction = "pca", dims = 1:20)
HCdat <- FindClusters(HCdat, resolution = 0.4)

table(HCdat$Condition)
Idents(HCdat) <- "Condition"
saveRDS(HCdat,file = "HC_Merged")









#============================= Merge RA HC (RA vs HC) ============================= 

HC <- readRDS("HC_Merged")
RA <- readRDS("RA_Merged")

obj.list<-list(HC,RA)
data.anchors<-FindIntegrationAnchors(obj.list, dims = 1:20)
M.dat<-IntegrateData(anchorset = data.anchors, dims = 1:20)
M.dat <- JoinLayers(M.dat,assay="RNA",layer="counts")
counts <- GetAssayData(M.dat, assay = "RNA",layer = "counts")

Tdat <- CreateSeuratObject(counts = counts,meta.data = M.dat@meta.data)
table(Tdat$Condition)

Tdat[["percent.mt"]]<- PercentageFeatureSet(Tdat, pattern = "^MT-")
Tdat <- subset(Tdat, subset= nFeature_RNA >200 & nFeature_RNA < 4000 & percent.mt <10)
Tdat <- NormalizeData(Tdat)
Tdat <- FindVariableFeatures(Tdat, selection.method = "vst", nFeature = 2000)
# Run PCA
all.genes <- rownames(Tdat)
Tdat <- ScaleData(Tdat, verbose = FALSE)
Tdat <- RunPCA(Tdat, npcs = 20, verbose = FALSE)
Tdat<- RunUMAP(Tdat, reduction = "pca", dims = 1:20)
Tdat <- FindNeighbors(Tdat, reduction = "pca", dims = 1:20)
Tdat <- FindClusters(Tdat, resolution = 0.4)

table(Tdat$Condition)
Idents(Tdat) <- "Condition"
saveRDS(Tdat,file = "RA_HC_Merged")

# Use find marker to get the difference matrix
b_markers <- FindMarkers(Tdat,
                         ident.1 = "RA",
                         ident.2 = "HC",
                         group.by ="Condition",
                         only.pos = FALSE)

b_markers
head(b_markers)
write.csv(x = b_markers, file = "123456_RA_VS_HC.csv") # Save it for volcano plot
