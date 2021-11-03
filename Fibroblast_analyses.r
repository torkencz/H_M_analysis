---
title: "R fib analysis"
output: html_notebook
---

Internal functions
```{r}
corrected_reads<-function(object,assay,regions,group_name,species=NULL,min=0,mid=0.1,max=0.2){

list_of_plots=list()
  
if(is.null(species)) {
species=c("mouse","human")
}
object_use=object
#lets look at corrected accessibility
  
require(ComplexHeatmap)  
require(scales)  
require(circlize)  
 
if(is.character(regions)){
str_test=regions
gr_test=StringToGRanges(regions)
}else
{
gr_test=regions
str_test=GRangesToString(regions)
}

for (sp in species){
#name of a factor in the object
print(sp)
object_sub<-subset(object_use,subset = species == paste0(sp))
groups<-as.factor(object_sub[[]][,group_name])
print("got here")

C<-GetAssayData(object_sub,slot = "counts",assay = assay)[str_test,]
C<-as.matrix((C > 0) + 0)
#aggregate
Cagg<-aggregate(t(C),by=list(groups),FUN = function(x){sum(x)/length(x)})
Cagg2<-as.data.frame(t(Cagg))
colnames(Cagg2)<-Cagg2[1,]
Cagg2<-Cagg2[-1,]
#turn everything into numeric
Cagg2[] <- lapply(Cagg2, function(x) {
   as.numeric(as.character(x))
})
#at this point we have the probability of accessibility for each cluster. To adjust for complexity we have to calculate the median number of sites accessible in each cluster we take the average of this and divide them by individual groups ()

median_nfeature_category<-aggregate(as.data.frame(object_sub$nFeature_peaks),by=list(groups),FUN=median)

median_nfeature_category$norm_complexity<-mean(log10(median_nfeature_category$`object_sub$nFeature_peaks`))/log10(median_nfeature_category$`object_sub$nFeature_peaks`)
norm_factor<-as.data.frame(median_nfeature_category$norm_complexity)
row.names(norm_factor)<-median_nfeature_category$Group.1
t(norm_factor)

temp<-apply(Cagg2,MARGIN = 1,FUN = function(x){data.frame(x*t(norm_factor))})
Cagg2<-ldply(temp, data.frame)
row.names(Cagg2)<-Cagg2$.id
Cagg2$.id<-NULL

plot_matrix<-as.matrix(Cagg2)

col_fun = colorRamp2(c(min, mid, max), c("blue", "white", "red"))
H<-Heatmap(plot_matrix,show_column_names = T,show_row_names = F,cluster_rows = T,col=col_fun)
list_of_plots[[paste0(sp)]]=H
}


return(list_of_plots)


}

```


```{r}
library(Seurat)
library(Signac)
integrated_H_M<-readRDS(file="/home/torkenczyk/H_M_chromatin_atlas/human_data/active_use_objects/Final_integration_all_bins_CCA_remove_removed.rds")

cellsF<-names(Idents(integrated_H_M)[Idents(integrated_H_M)=="Fibroblasts"])
cellsMF<-names(Idents(integrated_H_M)[Idents(integrated_H_M)=="Myofibroblasts"])

integrated_H_M<-subset(integrated_H_M,cells=c(cellsF,cellsMF))
integrated_H_M_fib <- RunUMAP(object = integrated_H_M, reduction = 'integrated_dr', dims = 2:100)


tissue_corr<-integrated_H_M_fib$tissue
tissue_corr[is.na(tissue_corr)]<-"unknown"
tissue_corr<-factor(tissue_corr)


levels(tissue_corr)<-c("Large intestine","Large intestine","Heart","Heart","Large intestine","Liver","Liver","Lung","Lung","Small intestine","Small intestine","unknown")

tissue_corr<-as.data.frame(tissue_corr)


integrated_H_M_fib<-AddMetaData(integrated_H_M_fib,tissue_corr)


DimPlot(integrated_H_M_fib,group.by = "tissue",split.by = "species",raster=F,label = TRUE)

DimPlot(integrated_H_M_fib,group.by = "tissue.replicate",split.by = "species",raster=F,label = TRUE)
DimPlot(integrated_H_M_fib,group.by = "ID1",split.by = "species",raster=F,label = TRUE)


library(cisTopic)

integrated_H_M_fib$celltype_corr<-Idents(integrated_H_M_fib)
rep_corr1<-integrated_H_M_fib$ID1[integrated_H_M_fib$species=="human"]
rep_corr1[is.na(rep_corr1)]<-"unknown"
rep_corr1<-as.data.frame(rep_corr1)
rep_corr2<-integrated_H_M_fib$tissue.replicate[integrated_H_M_fib$species=="mouse"]
rep_corr2<-as.data.frame(rep_corr2)
names(rep_corr1)="rep_corr"
names(rep_corr2)="rep_corr"
rep_corr<-rbind(rep_corr1,rep_corr2)
integrated_H_M_fib<-AddMetaData(integrated_H_M_fib,rep_corr)

saveRDS(integrated_H_M_fib,file="integrated_H_M_fib.rds")

load(file = "./active_use_objects/for_rahul/Shared.bw.species.DA.RData")
load(file = "./active_use_objects/for_rahul/Human_only.bw.species.DA.RData")
load(file = "./active_use_objects/for_rahul/Mouse_only.bw.species.DA.RData")

#create a combined cre object

human<-readRDS("./active_use_objects/Select_human_tissues.rds")
mouse<-readRDS("./active_use_objects/Select_mouse_tissues.rds")
 
human[["Deviation_DA_mouse"]]<-NULL
human[["Deviation_DA_human"]]<-NULL
human[["Deviation_DA_int"]]<-NULL
human[["chromvar"]]<-NULL
mouse[["Deviation_DA_mouse"]]<-NULL
mouse[["Deviation_DA_human"]]<-NULL
mouse[["Deviation_DA_int"]]<-NULL
mouse[["mm9"]]<-NULL
mouse[["mm10"]]<-NULL
mouse[["cre05"]]<-NULL
mouse[["cre99"]]<-NULL
mouse[["chromvar"]]<-NULL
human$species="human"
mouse$species="mouse"

H_M_merged_cre<-merge(human,mouse)

cellsF<-names(Idents(H_M_merged_cre)[Idents(H_M_merged_cre)=="Fibroblasts"])
cellsMF<-names(Idents(H_M_merged_cre)[Idents(H_M_merged_cre)=="Myofibroblasts"])
H_M_merged_cre_fib<-subset(H_M_merged_cre,cells=c(cellsF,cellsMF))

H_M_merged_cre_fib$celltype_corr<-Idents(H_M_merged_cre_fib)
rep_corr1<-H_M_merged_cre_fib$ID1[H_M_merged_cre_fib$species=="human"]
rep_corr1[is.na(rep_corr1)]<-"unknown"
rep_corr1<-as.data.frame(rep_corr1)
rep_corr2<-H_M_merged_cre_fib$tissue.replicate[H_M_merged_cre_fib$species=="mouse"]
rep_corr2<-as.data.frame(rep_corr2)
names(rep_corr1)="rep_corr"
names(rep_corr2)="rep_corr"
rep_corr<-rbind(rep_corr1,rep_corr2)
H_M_merged_cre_fib<-AddMetaData(H_M_merged_cre_fib,rep_corr)
tissue_corr<-H_M_merged_cre_fib$tissue
tissue_corr[is.na(tissue_corr)]<-"unknown"
tissue_corr<-factor(tissue_corr)
levels(tissue_corr)<-c("Large intestine","Large intestine","Heart","Heart","Large intestine","Liver","Liver","Lung","Lung","Small intestine","Small intestine","unknown")

tissue_corr<-as.data.frame(tissue_corr)


H_M_merged_cre_fib<-AddMetaData(H_M_merged_cre_fib,tissue_corr)



#now calculate a corrected 
sh<-corrected_reads(object = H_M_merged_cre_fib,assay = "cre",group_name = "tissue_corr",regions = shared_DA$Myofibroblasts$feature_hg38_universal,mid = 0.075,max = 0.15)

sh$human+sh$mouse


h_o<-corrected_reads(object = H_M_merged_cre_fib,assay = "cre",group_name = "tissue_corr",regions = human_DA$Myofibroblasts$feature_hg38_universal,mid = 0.075,max = 0.15)

h_o$human+h_o$mouse


m_o<-corrected_reads(object = H_M_merged_cre_fib,assay = "cre",group_name = "tissue_corr",regions = mouse_DA$Myofibroblasts$feature_hg38_universal,mid = 0.075,max = 0.15)
m_o$mouse+m_o$human


#do da where we quantify all fibrioblast specific sites


DimPlot(integrated_H_M_fib)


integrated_H_M_fib <- FindNeighbors(object = integrated_H_M_fib, reduction = 'integrated_dr', dims = 2:100)
integrated_H_M_fib <- FindClusters(object = integrated_H_M_fib, verbose = T, algorithm = 3,resolution = 0.1)

DimPlot(integrated_H_M_fib)

#from here

add<-H_M_merged_cre[["Row.names"]]
add$longer_name=row.names(add)
add$shorter_name=add$Row.names
add$Row.names<-NULL


clusters_add<-integrated_H_M_fib[["humanPeaks_snn_res.0.1"]][integrated_H_M$species=="mouse"]

clusters_add$shorter_name=row.names(clusters_add)


comb<-merge(clusters_add,add,by="shorter_name")

row.names(comb)<-comb$longer_name
comb$shorter_name<-NULL
comb$longer_name<-NULL

clusters_add<-integrated_H_M_fib$humanPeaks_snn_res.0.1[integrated_H_M_fib$species=="mouse"]
m_comb<-as.data.frame(clusters_add)
names(m_comb)<-"humanPeaks_snn_res.0.1"

annotation_renamed_rows<-rbind(comb,m_comb)



H_M_merged_cre_fib<-AddMetaData(H_M_merged_cre_fib,annotation_renamed_rows)

#now calculate a corrected 
sh<-corrected_reads(object = H_M_merged_cre_fib,assay = "cre",group_name = "humanPeaks_snn_res.0.1",regions = shared_DA$Myofibroblasts$feature_hg38_universal,mid = 0.075,max = 0.15)

sh$human+sh$mouse


h_o<-corrected_reads(object = H_M_merged_cre_fib,assay = "cre",group_name = "humanPeaks_snn_res.0.1",regions = human_DA$Myofibroblasts$feature_hg38_universal,mid = 0.1,max = 0.2)

h_o$human+h_o$mouse


m_o<-corrected_reads(object = H_M_merged_cre_fib,assay = "cre",group_name = "humanPeaks_snn_res.0.1",regions = mouse_DA$Myofibroblasts$feature_hg38_universal,mid = 0.075,max = 0.15)
m_o$mouse+m_o$human




#cistopic analysis
Fib_counts<-GetAssayData(integrated_H_M_fib,assay = "humanPeaks",slot = "counts")
row.names(Fib_counts)<-GRangesToString(StringToGRanges(row.names(Fib_counts)),sep = c(":","-"))
cisTopicObject <- createcisTopicObject(Fib_counts, project.name='Fib')
cellData_Fib<-integrated_H_M_fib[[c("celltype_corr","tissue_corr","species","rep_corr")]]
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = cellData_Fib)
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2:15, 20, 25, 35, 40, 45, 50), seed=123, nCores=10, addModels=FALSE)
# For WarpLDA
cisTopicObject@calc.params[['runWarpLDAModels']]$seed <- 123
cisTopicObject@calc.params[['runWarpLDAModels']]$iterations <- iterations
cisTopicObject@calc.params[['runWarpLDAModels']]$alpha <- alpha
cisTopicObject@calc.params[['runWarpLDAModels']]$alphaByTopic <- alphaByTopic
cisTopicObject@calc.params[['runWarpLDAModels']]$beta <- beta
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c("celltype_corr","tissue_corr","species","rep_corr"))
cellTopicHeatmap(cisTopicObject, method='Probability', colorBy=c("tissue_corr"))
cisTopicObject <- runUmap(cisTopicObject, target='cell')
par(mfrow=c(2,2))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c("celltype_corr","tissue_corr","species","rep_corr"), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)
par(mfrow=c(2,5))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)
cisTopicObject <- binarizedcisTopicsToCtx(cisTopicObject, genome='hg38')
saveRDS(cisTopicObject,file="Fib_cistopic.rds")
#


Fib_counts<-GetAssayData(integrated_H_M_fib,assay = "humanPeaks",slot = "counts")






DefaultAssay(Fib_ATAC_ac)<-"all"

enriched.motifs <- FindMotifs(
  object = Fib_ATAC_ac,
  features = GRangesToString(T6)
)
MotifPlot(
  object = Fib_ATAC_ac,
  motifs = head(rownames(enriched.motifs),n=20)
)

getBedFiles(cisTopicObject, path='output/cisTopics_asBed')







#remove NA cells and then do everything over again

integrated_H_M$Cell.Type[integrated_H_M$species]





```

