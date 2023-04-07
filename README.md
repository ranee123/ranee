# ranee
study
###学习irGSEA分析，注意如何把自己的细胞注释后的数据整理成类似pbmc3k的格式？？
#seurat_annotations如何自主注视到serat对象里面？？
library(Seurat)
library(SeuratData)
library(pbmc3k.SeuratData)##直接下载到本地再安装，网络太差了
# devtools::install_github('satijalab/seurat-data')
# view all available datasets
View(AvailableData())
# download 3k PBMCs from 10X Genomics
#InstallData("pbmc3k")##直接下载到本地再安装，网络太差了
# the details of pbmc3k.final
?pbmc3k.final

# loading dataset
data("pbmc3k.final")
pbmc3k.final <- UpdateSeuratObject(pbmc3k.final)
# plot
DimPlot(pbmc3k.final, reduction = "umap",
        group.by = "seurat_annotations",label = T) + NoLegend()
# set cluster to idents
Idents(pbmc3k.final) <- pbmc3k.final$seurat_annotations
str(Idents(pbmc3k.final) )
library(UCell)
library(irGSEA)
# calculate and integrate计算和整合
pbmc3k.final <- irGSEA.score(object = pbmc3k.final, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Homo sapiens", category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

Seurat::Assays(pbmc3k.final)
result.dge <- irGSEA.integrate(object = pbmc3k.final, 
                               group.by = "seurat_annotations",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))

class(result.dge)
#Visualization可视化

irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.heatmap.plot
