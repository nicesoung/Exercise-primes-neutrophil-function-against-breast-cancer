# Exercise-primes-neutrophil-function-against-breast-cancer
Analysis
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(plyr)
x.getPalette <- colorRampPalette(brewer.pal(12, "Paired"))





########
x.obj.filter1 <- readRDS("E:/nicesoung/Exercise_Neutrophil/GROUPSample/All_Blood/Harmony/R(5)0.3_Realfinal/Subtyping/Neutrophil/R(4)0.1_final/Seurat.subtype.RDS")
x.obj.filter2 <- readRDS("E:/nicesoung/Exercise_Neutrophil/GROUPSample/All_Bone/Harmony/R(7)0.1_Realfinal/Subtyping/Neutrophil/R(7)0.3_final/Seurat.subtype.RDS")
x.obj.filter3 <- readRDS("E:/nicesoung/Exercise_Neutrophil/GROUPSample/NP_Blood/Harmony/R(7)0.3_realfinal/Subtyping/Neutrophil/R(7)0.3_final/Seurat.subtype.RDS")
x.obj.filter4 <- readRDS("E:/nicesoung/Exercise_Neutrophil/GROUPSample/NP_Bone/Harmony/R(4)0.5_Realfinal/Subtyping/Neutrophil/R(5)0.1/Seurat.subtype.RDS")




########################### AddModuleScore
pathGeneTable = read.csv("E:/Research/Exercise/Scoring/N1_N2_gene.csv", header = T)
pathGeneTable


###################x.obj.filter1
for(i in 1:4){
pathGeneTable[i,pathGeneTable[i,] != ""]
l <- length(pathGeneTable[i,pathGeneTable[i,]!=""])
pathName <- pathGeneTable[i,1:2]
pathName <- pathName$V1
pathGeneSym <- as.character(pathGeneTable[i,3:l]) 
geneSet_list_1 <- list(pathGeneSym)
x.obj.filter1_AddModuleScore <- AddModuleScore(object = x.obj.filter1, features = c(geneSet_list_1))
Data <- data.frame(x.obj.filter1_AddModuleScore@meta.data$Cluster1)
rownames(Data) <- rownames(x.obj.filter1_AddModuleScore@meta.data)
colnames(Data) <- pathName
  if(i==1){
    Data1 <- Data
  }else{
    Data1 <- cbind(Data, Data1)
  }
}
head(Data)
head(Data1)
pathGeneSym
geneSet_list_1

x.obj.filter1_AddModuleScore@meta.data <-cbind(x.obj.filter1_AddModuleScore@meta.data,Data1)

saveRDS(x.obj.filter1_AddModuleScore, file=paste0('E:/Research/Exercise/Scoring/Seurat.AddModuleScore.Exercise_Neutrophil.All_Blood.Subtyping.Neutrophil_1.RDS'))




###################x.obj.filter2
for(i in 1:4){
  pathGeneTable[i,pathGeneTable[i,] != ""]
  l <- length(pathGeneTable[i,pathGeneTable[i,]!=""])
  pathName <- pathGeneTable[i,1:2]
  pathName <- pathName$V1
  pathGeneSym <- as.character(pathGeneTable[i,3:l]) 
  geneSet_list_1 <- list(pathGeneSym)
  x.obj.filter2_AddModuleScore <- AddModuleScore(object = x.obj.filter2, features = c(geneSet_list_1))
  Data <- data.frame(x.obj.filter2_AddModuleScore@meta.data$Cluster1)
  rownames(Data) <- rownames(x.obj.filter2_AddModuleScore@meta.data)
  colnames(Data) <- pathName
  if(i==1){
    Data1 <- Data
  }else{
    Data1 <- cbind(Data, Data1)
  }
}
head(Data)
head(Data1)
pathGeneSym
geneSet_list_1

x.obj.filter2_AddModuleScore@meta.data <-cbind(x.obj.filter2_AddModuleScore@meta.data,Data1)

saveRDS(x.obj.filter2_AddModuleScore, file=paste0('E:/Research/Exercise/Scoring/Seurat.AddModuleScore.Exercise_Neutrophil.All_Bone.Subtyping.Neutrophil_1.RDS'))





###################x.obj.filter3
for(i in 1:4){
  pathGeneTable[i,pathGeneTable[i,] != ""]
  l <- length(pathGeneTable[i,pathGeneTable[i,]!=""])
  pathName <- pathGeneTable[i,1:2]
  pathName <- pathName$V1
  pathGeneSym <- as.character(pathGeneTable[i,3:l]) 
  geneSet_list_1 <- list(pathGeneSym)
  x.obj.filter3_AddModuleScore <- AddModuleScore(object = x.obj.filter3, features = c(geneSet_list_1))
  Data <- data.frame(x.obj.filter3_AddModuleScore@meta.data$Cluster1)
  rownames(Data) <- rownames(x.obj.filter3_AddModuleScore@meta.data)
  colnames(Data) <- pathName
  if(i==1){
    Data1 <- Data
  }else{
    Data1 <- cbind(Data, Data1)
  }
}
head(Data)
head(Data1)
pathGeneSym
geneSet_list_1

x.obj.filter3_AddModuleScore@meta.data <-cbind(x.obj.filter3_AddModuleScore@meta.data,Data1)

saveRDS(x.obj.filter3_AddModuleScore, file=paste0('E:/Research/Exercise/Scoring/Seurat.AddModuleScore.Exercise_Neutrophil.NP_Blood.Subtyping.Neutrophil_1.RDS'))





###################x.obj.filter4
for(i in 1:4){
  pathGeneTable[i,pathGeneTable[i,] != ""]
  l <- length(pathGeneTable[i,pathGeneTable[i,]!=""])
  pathName <- pathGeneTable[i,1:2]
  pathName <- pathName$V1
  pathGeneSym <- as.character(pathGeneTable[i,3:l]) 
  geneSet_list_1 <- list(pathGeneSym)
  x.obj.filter4_AddModuleScore <- AddModuleScore(object = x.obj.filter4, features = c(geneSet_list_1))
  Data <- data.frame(x.obj.filter4_AddModuleScore@meta.data$Cluster1)
  rownames(Data) <- rownames(x.obj.filter4_AddModuleScore@meta.data)
  colnames(Data) <- pathName
  if(i==1){
    Data1 <- Data
  }else{
    Data1 <- cbind(Data, Data1)
  }
}
head(Data)
head(Data1)
pathGeneSym
geneSet_list_1

x.obj.filter4_AddModuleScore@meta.data <-cbind(x.obj.filter4_AddModuleScore@meta.data,Data1)

saveRDS(x.obj.filter4_AddModuleScore, file=paste0('E:/Research/Exercise/Scoring/Seurat.AddModuleScore.Exercise_Neutrophil.NP_Bone.Subtyping.Neutrophil_1.RDS'))






All_Blood <- x.obj.filter1_AddModuleScore
All_Bone <- x.obj.filter2_AddModuleScore
NP_Blood <- x.obj.filter3_AddModuleScore
NP_Bone <- x.obj.filter4_AddModuleScore



All_Blood@meta.data
write.csv(All_Blood@meta.data,"E:/Research/Exercise/Scoring/All_Blood_raw.csv")
write.csv(All_Bone@meta.data,"E:/Research/Exercise/Scoring/All_Bone_raw.csv")
write.csv(NP_Blood@meta.data,"E:/Research/Exercise/Scoring/NP_Blood_raw.csv")
write.csv(NP_Bone@meta.data,"E:/Research/Exercise/Scoring/NP_Bone_raw.csv")



##########ALL_Blood
VlnPlot(object= All_Blood, c('N1_Anti_tumor_Markers_for_Blood_derived_Neutrophils'), group.by = 'Batch')
png(filename = paste("E:/Research/Exercise/Scoring/PNG/All_Blood_N1_vln.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(VlnPlot(object= All_Blood, c('N1_Anti_tumor_Markers_for_Blood_derived_Neutrophils'), group.by = 'Batch'))
dev.off()

VlnPlot(object= All_Blood, c('N2_Pro_tumor_Markers_for_Blood_derived_Neutrophils'), group.by = 'Batch')
png(filename = paste("E:/Research/Exercise/Scoring/PNG/All_Blood_N2_vln.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(VlnPlot(object= All_Blood, c('N2_Pro_tumor_Markers_for_Blood_derived_Neutrophils'), group.by = 'Batch'))
dev.off()

FeaturePlot(All_Blood, features =  c('N1_Anti_tumor_Markers_for_Blood_derived_Neutrophils'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
png(filename = paste("E:/Research/Exercise/Scoring/PNG/All_Blood_N1_Featureplot.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(FeaturePlot(All_Blood, features =  c('N1_Anti_tumor_Markers_for_Blood_derived_Neutrophils'))& 
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
dev.off()

FeaturePlot(All_Blood, features =  c('N2_Pro_tumor_Markers_for_Blood_derived_Neutrophils'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
png(filename = paste("E:/Research/Exercise/Scoring/PNG/All_Blood_N2_Featureplot.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(FeaturePlot(All_Blood, features =  c('N2_Pro_tumor_Markers_for_Blood_derived_Neutrophils'))& 
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
dev.off()



###All_BLood N1 marker
FeaturePlot(All_Blood, features =  c('Nos2'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(All_Blood, features =  c('Tnf'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(All_Blood, features =  c('Cd86'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(All_Blood, features =  c('Cxcl9'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(All_Blood, features =  c('Cxcl10'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(All_Blood, features =  c('Itgam'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(All_Blood, features =  c('Icam1'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


###All_BLood N2 marker
FeaturePlot(All_Blood, features =  c('Arg1'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(All_Blood, features =  c('Cxcr4'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(All_Blood, features =  c('Il10'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(All_Blood, features =  c('Tgfb1'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(All_Blood, features =  c('S100a8'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
FeaturePlot(All_Blood, features =  c('S100a9'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


VlnPlot(object= All_Blood, c('S100a8'), group.by = 'Batch')
VlnPlot(object= All_Blood, c('S100a9'), group.by = 'Batch')






##########NP_Blood
VlnPlot(object= NP_Blood, c('N1_Anti_tumor_Markers_for_Blood_derived_Neutrophils'), group.by = 'Batch')
png(filename = paste("E:/Research/Exercise/Scoring/PNG/NP_Blood_N1_vln.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(VlnPlot(object= NP_Blood, c('N1_Anti_tumor_Markers_for_Blood_derived_Neutrophils'), group.by = 'Batch'))
dev.off()

VlnPlot(object= NP_Blood, c('N2_Pro_tumor_Markers_for_Blood_derived_Neutrophils'), group.by = 'Batch')
png(filename = paste("E:/Research/Exercise/Scoring/PNG/NP_Blood_N2_vln.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(VlnPlot(object= NP_Blood, c('N2_Pro_tumor_Markers_for_Blood_derived_Neutrophils'), group.by = 'Batch'))
dev.off()

FeaturePlot(NP_Blood, features =  c('N1_Anti_tumor_Markers_for_Blood_derived_Neutrophils'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
png(filename = paste("E:/Research/Exercise/Scoring/PNG/NP_Blood_N1_Featureplot.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(FeaturePlot(NP_Blood, features =  c('N1_Anti_tumor_Markers_for_Blood_derived_Neutrophils'))& 
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
dev.off()

FeaturePlot(NP_Blood, features =  c('N2_Pro_tumor_Markers_for_Blood_derived_Neutrophils'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
png(filename = paste("E:/Research/Exercise/Scoring/PNG/NP_Blood_N2_Featureplot.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(FeaturePlot(NP_Blood, features =  c('N2_Pro_tumor_Markers_for_Blood_derived_Neutrophils'))& 
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
dev.off()





##########All_Bone
VlnPlot(object= All_Bone, c('N1_Anti_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'), group.by = 'Batch')
png(filename = paste("E:/Research/Exercise/Scoring/PNG/All_Bone_N1_vln.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(VlnPlot(object= All_Bone, c('N1_Anti_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'), group.by = 'Batch'))
dev.off()

VlnPlot(object= All_Bone, c('N2_Pro_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'), group.by = 'Batch')
png(filename = paste("E:/Research/Exercise/Scoring/PNG/All_Bone_N2_vln.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(VlnPlot(object= All_Bone, c('N2_Pro_tumor_Markers_for_Blood_derived_Neutrophils'), group.by = 'Batch'))
dev.off()

FeaturePlot(All_Bone, features =  c('N1_Anti_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
png(filename = paste("E:/Research/Exercise/Scoring/PNG/All_Bone_N1_Featureplot.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(FeaturePlot(All_Bone, features =  c('N1_Anti_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'))& 
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
dev.off()

FeaturePlot(All_Bone, features =  c('N2_Pro_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
png(filename = paste("E:/Research/Exercise/Scoring/PNG/All_Bone_N2_Featureplot.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(FeaturePlot(All_Bone, features =  c('N2_Pro_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'))& 
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
dev.off()


##########NP_Bone
VlnPlot(object= NP_Bone, c('N1_Anti_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'), group.by = 'Batch')
png(filename = paste("E:/Research/Exercise/Scoring/PNG/NP_Bone_N1_vln.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(VlnPlot(object= NP_Bone, c('N1_Anti_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'), group.by = 'Batch'))
dev.off()

VlnPlot(object= NP_Bone, c('N2_Pro_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'), group.by = 'Batch')
png(filename = paste("E:/Research/Exercise/Scoring/PNG/NP_Bone_N2_vln.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(VlnPlot(object= NP_Bone, c('N2_Pro_tumor_Markers_for_Blood_derived_Neutrophils'), group.by = 'Batch'))
dev.off()

FeaturePlot(NP_Bone, features =  c('N1_Anti_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
png(filename = paste("E:/Research/Exercise/Scoring/PNG/NP_Bone_N1_Featureplot.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(FeaturePlot(NP_Bone, features =  c('N1_Anti_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'))& 
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
dev.off()

FeaturePlot(NP_Bone, features =  c('N2_Pro_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
png(filename = paste("E:/Research/Exercise/Scoring/PNG/NP_Bone_N2_Featureplot.png", sep=""), width=2000, height=2000,res=200, type="cairo")
print(FeaturePlot(NP_Bone, features =  c('N2_Pro_tumor_Markers_for_Bone_Marrow_derived_Neutrophils'))& 
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
dev.off()
