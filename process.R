
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gc()
memory.limit(9999999999)
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)

setwd("D:/analysis/forpublication/baso_merge_new")
set.seed(123)
getwd()
library(Seurat)#导入R包
data_dir <- "D:/analysis/forpublication/baso/filtered_feature_bc_matrix"
data <- Read10X(data.dir = data_dir,)#设置文件存储的位置，然后读取
baso1 <- CreateSeuratObject(counts = data,project = "baso1")

data_dir <- "D:/analysis/forpublication/baso2/filtered_feature_bc_matrix"
data <- Read10X(data.dir = data_dir,)#设置文件存储的位置，然后读取
baso2 <- CreateSeuratObject(counts = data,project = "baso2")

YDL<-merge(baso1,baso2)

pdf("质控.pdf")
YDL[["percent.mt"]] <- PercentageFeatureSet(object = YDL, pattern = "^mt-")
###展示基因及线粒体百分比（这里将其进行标记并统计其分布频率，"nFeature_RNA"为基因数，"nCount_RNA"为UMI数，"percent.mt"为线粒体占比）
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = YDL, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

##过滤细胞：根据上面小提琴图中基因数"nFeature_RNA"和线粒体数"percent.mt"，分别设置过滤参数，这里基因数 200-4000，线粒体百分比为小于 5%，保留gene数大于200小于2500的细胞；目的是去掉空GEMs和1个GEMs包含2个以上细胞的数据；而保留线粒体基因的转录本数低于5%的细胞,为了过滤掉死细胞等低质量的细胞数据。
YDL <- subset(x = YDL, subset = nFeature_RNA > 200 &nFeature_RNA < 4000 &nCount_RNA<40000 &nCount_RNA>1000 & percent.mt < 15)    #对数据进行过滤
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = YDL, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
YDL
table(Idents(YDL))
summary(YDL@meta.data$nCount_RNA)
summary(YDL@meta.data$nFeature_RNA)
summary(YDL@meta.data$percent.mt)
dev.off()


#对数据进行标准化
##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000
YDL <- NormalizeData(object = YDL, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
YDL <- FindVariableFeatures(object = YDL, selection.method = "vst", nfeatures = 2000)
YDL <- ScaleData(object = YDL, features = rownames(YDL))
#线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
YDL=RunPCA(object= YDL,npcs = 20,pc.genes=VariableFeatures(object = YDL))     #PCA分析
ElbowPlot(YDL)#选择top20个PC
pcSelect=10
YDL <- FindNeighbors(object = YDL, dims = 1:pcSelect)                #计算邻接距离
##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
YDL <- FindClusters(object = YDL, resolution = 0.2)                  #对细胞分组,优化标准模块化
##使用Idents（）函数可查看不同细胞的分群；
head(Idents(YDL), 5)

##Seurat提供了几种非线性降维的方法进行数据可视化（在低维空间把相似的细胞聚在一起）,比如UMAP和t-SNE,运行UMAP需要先安装'umap-learn'包，这里不做介绍，两种方法都可以使用，但不要混用，如果混用，后面的结算结果会将先前的聚类覆盖掉，只能保留一个。

##这里采用基于TSNE的聚类方法。
YDL <- RunTSNE(object = YDL, dims = 1:pcSelect,check_duplicates = FALSE)                      #TSNE聚类
pdf(file="TSNE_baso.pdf",width=6.5,height=6)
TSNEPlot(object = YDL, pt.size = 2, label = TRUE)    #TSNE可视化

#另一个可视化的方法
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
dev.off()
write.table(YDL$seurat_clusters,file="tsneCluster_baso.txt",quote=F,sep="\t",col.names=F)


##这里采用基于umap的聚类方法
#这里采用基于图论的聚类方法
YDL <- RunUMAP(object = YDL, dims = 1:pcSelect)                      #umap聚类
pdf(file="umap__baso.pdf",width=6.5,height=6)
UMAPPlot(object = YDL, pt.size = 1.5, label = TRUE)    #umap可视化
#另一个可视化的方法
DimPlot(object=YDL,label = TRUE,reduction="umap")
##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；
dev.off()

##细胞周期归类

YDL<- CellCycleScoring(object = YDL, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

head(x = YDL@meta.data)
pdf(file="CellCycle_baso.pdf",width=6.5,height=6)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)
dev.off()


library(rSuperCT)
library(ggplot2)
pred_obj <- ImportData(YDL)
dir.create('./models', showWarnings = FALSE)
pred_obj <- PredCellTypes(pred_obj, species = 'mouse', model = 'generic_38celltypes',results.dir = './models')
table(pred_obj@meta.data$pred_types)# 举例来展示不同细胞的比例
g <- plotHist(pred_obj) + scale_fill_manual(values = rep('blue', 13))
write.csv(g$data, 'pred_types.hist.csv',row.names = F)
#选择特定类型的细胞
Idents(YDL) <- pred_obj@meta.data$pred_types
#YDL <- subset(YDL, idents = 'Redblood')
#save(YDL, file = 'YDL.Redblood.Rdata')
Idents(YDL) <- YDL@meta.data$seurat_clusters
YDL <- subset(YDL, idents = c("0","1","2","3","4","5"))
#另一个可视化的方法
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)



set.seed(123)

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

library(DoubletFinder)

library(Seurat)

library(dplyr)

library(Matrix)

library(methods)

library(RColorBrewer)

YDL#读入处理好的seurat对象
sweep.res.list <- paramSweep_v3(YDL, PCs = 1:20, sct = FALSE)
head(sweep.res.list)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
head(sweep.res.list)
bcmvn <- find.pK(sweep.stats)

annotations <- YDL@meta.data$seurat_clusters

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------

homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- YDL@meta.data$ClusteringResults

nExp_poi <- round(0.075*nrow(YDL@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))



## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------

YDL <- doubletFinder_v3(YDL, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(YDL@meta.data)
YDL <- doubletFinder_v3(YDL, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_398", sct = FALSE)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "orig.ident")    #TSNE可视化
YDL_NEW <- subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$DF.classifications_0.25_0.09_398=="Singlet",]))
DimPlot(YDL_NEW,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "orig.ident")
DimPlot(YDL_NEW,reduction = "tsne",label = TRUE,pt.size = 1.5)

YDL_NEW_1 <- subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$DF.classifications_0.25_0.09_398=="Doublet",]))
DimPlot(YDL_NEW_1,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "orig.ident")
DimPlot(YDL_NEW_1,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by="DF.classifications_0.25_0.09_398")



saveRDS(YDL_NEW,"baso_merge_SINGLET.RDS")

pdf("DoubletFinder处理前后.pdf")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by="DF.classifications_0.25_0.09_398")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "orig.ident")    #TSNE可视化
dev.off()






Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gc()
memory.limit(9999999999)
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)

setwd("D:/analysis/forpublication/baso_merge_new")
YDL<-readRDS("baso_merge_SINGLET.RDS")

pdf(file="细胞周期_baso.pdf",width=6.5,height=6)
TSNEPlot(object = YDL, pt.size = 2, label = TRUE)    #TSNE可视化

#另一个可视化的方法
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "orig.ident")

DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,group.by = "orig.ident")


DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)
dev.off()

pdf("umi.pdf")
#绘制UMI的分布图
library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = log(nCount_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1
p2
plot_grid(p1,p2)
library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nCount_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5)

DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)


#绘制基因的分布图
library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","nFeature_RNA"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = log(nFeature_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1
p2
plot_grid(p1,p2)
library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nFeature_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nFeature_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5)

DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)



#绘制线粒体的分布图
library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","percent.mt"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = percent.mt))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1
p2
plot_grid(p1,p2)
library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","percent.mt"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = percent.mt))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5)

DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)

dev.off()


library("Seurat")
library("ggplot2")
gg <- TSNEPlot(YDL)
col<-ggplot_build(gg)$data
col<-as.data.frame(col)
table(col$colour)
table(YDL$seurat_clusters)

pdf("不同样本的cluster比例_小数点.pdf")
# Color palette
#colors <- c("#F8766D","#E68613","#CD9600","#ABA300","#7CAE00","#0CB702","#00BE67","#00C19A","#00BFC4","#00B8E7","#00A9FF","#8494FF","#C77CFF","#ED68ED","#FF61CC","#FF68A1")
colors <- c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")

DimPlot(YDL,cols =colors,reduction="tsne")
a<-as.data.frame(table(YDL@meta.data$seurat_clusters))
a<-a[c(1:6),]
# 数据准备
info = a$Freq
# 命名
names = as.character(a$Var1)
#names = c("0","1","2","3","4","5","6","7")
# 涂色（可选）
#cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC","#3f90CC","#22B17C","#FFC93E")
# 计算百分比
piepercent = paste(round(100*info/sum(info),2), "%")
# 绘图
pie(info, labels=piepercent, main = "total cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.8, fill=colors)
dev.off()


pdf("cell cycle ratio.pdf")

DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")

a<-as.data.frame(table(YDL@meta.data$Phase))
# 数据准备
info = a$Freq
# 命名
names = c("G1","G2M","S")
# 涂色（可选）
cols = c("#F8766D","#0CB702","#00A9FF")
# 计算百分比
piepercent = paste(round(100*info/sum(info)), "%")
# 绘图
pie(info, labels=piepercent, main = "cell cycle ratio", col=cols, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.8, fill=cols)
dev.off()


my_levels <- c(3,2,1,0,4,5)

factor(Idents(YDL), levels= my_levels) 
Idents(YDL) <- factor(Idents(YDL), levels= my_levels)
DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5,cols = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3"))

Idents(YDL) <- YDL@meta.data$seurat_clusters

YDL.AVERAGE<-AverageExpression(object =YDL,return.seurat=F)
YDL.AVERAGE<-as.data.frame(YDL.AVERAGE)

colnames(YDL.AVERAGE)<-c("cluster3","cluster2","cluster1","cluster0","cluster4","cluster5")

Idents(YDL) <- YDL@meta.data$orig.ident

YDL.AVERAGE<-AverageExpression(object =YDL,return.seurat=F)
YDL.AVERAGE<-as.data.frame(YDL.AVERAGE)

colnames(YDL.AVERAGE)<-c("baso1","baso2")




pdf("pearson.pdf")
results<-YDL.AVERAGE
cov(results)
pheatmap::pheatmap(cov(results))#,filename = "heatmap.pdf")
cor(results,method = "pearson")
pheatmap::pheatmap(cor(results,method = "pearson"),
                   cluster_rows = T,display_numbers = T,
                   cluster_cols = T)#,,filename = "heatmap.pdf")
pheatmap::pheatmap(cor(results,method = "spearman"),
                   cluster_rows = T,display_numbers = T,
                   cluster_cols = T)#,,filename = "heatmap.pdf")
pheatmap::pheatmap(cor(results,method = "spearman"),clustering_method = "single",
                   cluster_rows = T,display_numbers = T,
                   cluster_cols = T)#,,filename = "heatmap.pdf")
pheatmap::pheatmap(cor(results,method = "spearman"),clustering_method = "single",
                   cluster_rows = F,display_numbers = T,
                   cluster_cols = F)#,,filename = "heatmap.pdf")

pheatmap::pheatmap(cor(results,method = "spearman"),clustering_method = "ward.D2",
                   cluster_rows = T,display_numbers = T,
                   cluster_cols = T,color = colorRampPalette(colors = c("white","red"))(100))#,,filename = "heatmap.pdf")


dev.off()



YDL.markers <- FindAllMarkers(YDL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
#YDL.markers <- FindAllMarkers(YDL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 't')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run evYDL time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
YDL.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
##存储marker

write.csv(YDL.markers,file="allmarker_baso_321045.csv")
YDL.markers<-read.csv("allmarker_baso_321045.csv")
#绘制分cluster的热图
top10 <- YDL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#绘制marker在各个cluster的热图
pdf(file="tsneHeatmap_all.pdf",width=12,height=9)
DoHeatmap(object = YDL, features = top10$gene) + NoLegend()
DoHeatmap(subset(YDL, downsample = 100), features = top10$gene, size = 3)+ NoLegend()
DoHeatmap(object = YDL, features = YDL.markers$gene) + NoLegend()
DoHeatmap(subset(YDL, downsample = 100), features = YDL.markers$gene, size = 3)+ NoLegend()
dev.off()

dim(YDL.markers)
table(YDL.markers$cluster)

my_levels <- c(3,2,1,0,4,5)

factor(Idents(YDL), levels= my_levels) 
Idents(YDL) <- factor(Idents(YDL), levels= my_levels)
DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5,cols = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3"))
DoHeatmap(object = YDL, features = top10$gene,group.colors = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) + NoLegend()
DoHeatmap(object = YDL, features = YDL.markers$gene,group.colors = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) + NoLegend()
DoHeatmap(object = subset(YDL, downsample = 100), features = YDL.markers$gene,group.colors = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) + NoLegend()
DoHeatmap(object = YDL, features = top10$gene,group.colors = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) + NoLegend()
DoHeatmap(object = YDL, features = YDL.markers$gene,group.colors = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) + NoLegend()
DoHeatmap(object = subset(YDL, downsample = 100), features = YDL.markers$gene,group.colors = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) + NoLegend()



gene<- YDL.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

DotPlot(object = YDL, features = unique(gene$gene))+ RotatedAxis()
DotPlot(object = YDL, features =unique(gene$gene))+ RotatedAxis()
DotPlot(object = YDL, features =unique(gene$gene))+ RotatedAxis()+ coord_flip()

#加边框
DotPlot(object = YDL, features =unique(gene$gene))+ coord_flip()+ RotatedAxis()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)



YDL.AVERAGE<-AverageExpression(object = YDL,return.seurat=F)
YDL.AVERAGE<-as.data.frame(YDL.AVERAGE)

write.csv(YDL.AVERAGE,file = "AVERAGE_baso.csv")
data<-YDL.AVERAGE

merge_tf<-top10$gene
#挑选部分感兴趣
my.regulons <- merge_tf
#删掉所有列上都重复的
newdata<-data[c(my.regulons),]
newdata<-na.omit(newdata)
colnames(newdata)<-c("baso_cluster3","baso_cluster2","baso_cluster1","baso_cluster0","baso_cluster4","baso_cluster5")
#低值为蓝色，高值为红色，中间值为白色：
#pdf("fig1g.pdf")
pheatmap(newdata,fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))

pheatmap(newdata,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))




# 加载包
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))




pheatmap(newdata,scale = "row",fontsize = 6,filename = "new2.pdf",width = 10,height = 100,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()


library(monocle)
#准备monocle分析需要的文件
monocle.matrix=as.matrix(YDL@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="monocleMatrix.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(YDL@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="monocleGene.txt",quote=F,sep="\t",row.names=F)
write.table(YDL.markers,file="monocleMarkers.txt",sep="\t",row.names=F,quote=F)


#设置工作目录
monocle.matrix=read.table("monocleMatrix.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.sample=read.table("monocleSample.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("monocleGene.txt",sep="\t",header=T,row.names=1,check.names=F)
marker=read.table("monocleMarkers.txt",sep="\t",header=T,check.names=F)

#将Seurat结果转换为monocle需要的细胞矩阵，细胞注释表和基因注释表表
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)


#给其中一列数据重命名
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])
# saveRDS(cds,"cds.rds")
# rm(list = ls())  
# marker<-read.csv("allmarker.csv")
# cds<-readRDS("cds.rds")
#伪时间分析流程
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, marker$gene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')
cds <- orderCells(cds,reverse = T)
pdf(file="cluster.trajectory_SINGLET.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
plot_cell_trajectory(cds,color_by="Cluster")+facet_wrap(~Cluster,nrow=3,ncol = 2)
plot_cell_trajectory(cds,color_by="orig.ident")+facet_wrap(~orig.ident,nrow=2,ncol = 5)
plot_cell_trajectory(cds,color_by="Cluster")+facet_wrap(~orig.ident,nrow=2,ncol = 3)

plot_cell_trajectory(cds, color_by = "YDL@active.ident")
plot_cell_trajectory(cds, color_by="Pseudotime", show_backbone=FALSE)
# 可以很明显看到细胞的发育轨迹 
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State") +facet_wrap(~State, nrow = 1)
dev.off()

saveRDS(cds,"cds_SINGLET.rds")





detach("package:monocle")
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)


#YDL<-readRDS("SINGLET.RDS")

pdf("momocle3结果.pdf")
DimPlot(YDL, reduction = "umap",pt.size = 1.5,label = T)
###Pseudotime monocle3
cds <- as.cell_data_set(YDL)
cds <- cluster_cells(cds)
head(pData(cds))
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)


stem2<-rownames(pData(cds)[which(pData(cds)$ident %in% c('3')),])
cds <- order_cells(cds, root_cells = stem2)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
dev.off()







library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)
#合并四个时期的BP，以热图展示
a <- read.csv("allmarker_baso_321045.csv")
b<-a[a$cluster=="0","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene_cluster0 <- eg[,2]
gene_cluster0<-as.data.frame(gene_cluster0)
colnames(gene_cluster0)<-c("cluster0")
head(gene_cluster0)
gene_cluster0<-c(gene_cluster0)


a <- read.csv("allmarker_baso_321045.csv")
b<-a[a$cluster=="1","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene_cluster1 <- eg[,2]
gene_cluster1<-as.data.frame(gene_cluster1)
colnames(gene_cluster1)<-c("cluster1")
head(gene_cluster1)
gene_cluster1<-c(gene_cluster1)

a <- read.csv("allmarker_baso_321045.csv")
b<-a[a$cluster=="2","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene_cluster2 <- eg[,2]
gene_cluster2<-as.data.frame(gene_cluster2)
colnames(gene_cluster2)<-c("cluster2")
head(gene_cluster2)
gene_cluster2<-c(gene_cluster2)

a <- read.csv("allmarker_baso_321045.csv")
b<-a[a$cluster=="3","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene_cluster3 <- eg[,2]
gene_cluster3<-as.data.frame(gene_cluster3)
colnames(gene_cluster3)<-c("cluster3")
head(gene_cluster3)
class(gene_cluster3)
gene_cluster3<-c(gene_cluster3)

a <- read.csv("allmarker_baso_321045.csv")
b<-a[a$cluster=="4","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene_cluster4 <- eg[,2]
gene_cluster4<-as.data.frame(gene_cluster4)
colnames(gene_cluster4)<-c("cluster4")
head(gene_cluster4)
class(gene_cluster4)
gene_cluster4<-c(gene_cluster4)

a <- read.csv("allmarker_baso_321045.csv")
b<-a[a$cluster=="5","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene_cluster5 <- eg[,2]
gene_cluster5<-as.data.frame(gene_cluster5)
colnames(gene_cluster5)<-c("cluster5")
head(gene_cluster5)
class(gene_cluster5)
gene_cluster5<-c(gene_cluster5)


#data<-structure(list(cluster3 =gene_cluster3,cluster2 =gene_cluster2,cluster1 =gene_cluster1,cluster0 =gene_cluster0, cluster4=gene_cluster4,  cluster5=gene_cluster5))
data<-list(cluster3 =gene_cluster3$cluster3,cluster2 =gene_cluster2$cluster2,cluster1 =gene_cluster1$cluster1,cluster0 =gene_cluster0$cluster0, cluster4=gene_cluster4$cluster4,  cluster5=gene_cluster5$cluster5)

lapply(data, head)


head(as.data.frame(ck))

ck <- compareCluster(geneCluster = data,OrgDb = org.Mm.eg.db, fun = "enrichGO", pvalueCutoff=0.05)
write.csv(ck,"go.csv")


dim(ck)
dotplot(ck, showCategory =50)
ego2 <- simplify(ck,cutoff=0.7,by="p.adjust",select_fun=min)  #去除冗余，可以调整cutoff值
dim(ego2)
dotplot(ck, showCategory =50)
dotplot(ego2, showCategory =50)
dotplot(ego2, showCategory =15)

ck_enrichKEGG <- compareCluster(geneCluster = data, fun = "enrichKEGG",organism="mmu")
dotplot(ck_enrichKEGG, showCategory =15)
write.csv(ck_enrichKEGG,"enrichKEGG.csv")

#visualize the result using dotplot method.

pdf("多样本富集分析-1.pdf",width=10,height=12)
dotplot(ck, showCategory =50)
dotplot(ck, showCategory =30,split="ONTOLOGY") 
dev.off()


pdf("基因表达水平.pdf",width=25,height=12)
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Hbb-bt","Hbb-bs","Bpgm", "Alas2","Fam220a",
                                                                         "Tcp11l2","Cd36","Xpo7", "Tbcel","Tent5c",
                                                                        "Slc4a1", "Gpx1", "Slc25a37", "Rbm38", "Ncoa4"),cols = c("gray", "red"),ncol = 5)#actin

dev.off()


pdf("细胞周期比例.pdf")
a<-as.data.frame(table(YDL@meta.data$Phase))
# 数据准备
info = a$Freq
# 命名
names = c("G1","G2M","S")
# 涂色（可选）
cols = c("#F8766D","#0CB702","#00A9FF")
# 计算百分比
piepercent = paste(round(100*info/sum(info),2), "%")
# 绘图
pie(info, labels=piepercent, main = "cell cycle ratio", col=cols, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.8, fill=cols)


DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
dev.off()




dir.create("./pro")
setwd("./pro")
YDL<-readRDS("D:/analysis/forpublication/baso_merge_new/baso_merge_SINGLET.RDS")
wt_ctl_bm1<-readRDS("D:/analysis/forpublication/summary/pro_own/pro_SINGLET.RDS")
wt_ctl_bm1@meta.data$orig.ident<-c("pro")
YDL@meta.data$orig.ident<-c("baso")

head(YDL@meta.data)

YDL<-merge(YDL,wt_ctl_bm1)
Idents(YDL) <- YDL@meta.data$orig.ident
##计算每个细胞的线粒体基因转录本数的百分比（%）,使用[[ ]] 操作符存放到metadata中，mit-开头的为线粒体基因
YDL[["percent.mt"]] <- PercentageFeatureSet(object = YDL, pattern = "^mt-")
###展示基因及线粒体百分比（这里将其进行标记并统计其分布频率，"nFeature_RNA"为基因数，"nCount_RNA"为UMI数，"percent.mt"为线粒体占比）
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


head(YDL@meta.data)
summary(YDL@meta.data$nCount_RNA)
summary(YDL@meta.data$nFeature_RNA)
summary(YDL@meta.data$percent.mt)


##过滤细胞：根据上面小提琴图中基因数"nFeature_RNA"和线粒体数"percent.mt"，分别设置过滤参数，这里基因数 200-4000，线粒体百分比为小于 5%，保留gene数大于200小于2500的细胞；目的是去掉空GEMs和1个GEMs包含2个以上细胞的数据；而保留线粒体基因的转录本数低于5%的细胞,为了过滤掉死细胞等低质量的细胞数据。
# YDL <- subset(x = YDL, subset = nFeature_RNA > 200 &nFeature_RNA < 4000 &nCount_RNA<40000 &nCount_RNA>1000 & percent.mt < 15)    #对数据进行过滤
# VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(Idents(YDL))
summary(YDL@meta.data$nCount_RNA)
summary(YDL@meta.data$nFeature_RNA)
summary(YDL@meta.data$percent.mt)
# > summary(YDL@meta.data$nCount_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1062   17655   20865   21095   24848   39922 
# > summary(YDL@meta.data$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 201    1575    1962    1867    2298    3985 
# > summary(YDL@meta.data$percent.mt)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.03366  0.37025  0.58095  0.95798  0.99198 13.64552 


# YDL<-readRDS("pro_baso_SINGLET.RDS")
# saveRDS(YDL,"pro_baso_SINGLET.RDS")

#对数据进行标准化
##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000
YDL <- NormalizeData(object = YDL, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
YDL <- FindVariableFeatures(object = YDL, selection.method = "vst", nfeatures = 2000)
YDL <- ScaleData(object = YDL, features = rownames(YDL))
#线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
YDL=RunPCA(object= YDL,npcs = 20,pc.genes=VariableFeatures(object = YDL))     #PCA分析
ElbowPlot(YDL)#选择top20个PC
pcSelect=20
YDL <- FindNeighbors(object = YDL, dims = 1:pcSelect)                #计算邻接距离
##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
YDL <- FindClusters(object = YDL, resolution = 0.3)                  #对细胞分组,优化标准模块化
##使用Idents（）函数可查看不同细胞的分群；
head(Idents(YDL), 5)

##Seurat提供了几种非线性降维的方法进行数据可视化（在低维空间把相似的细胞聚在一起）,比如UMAP和t-SNE,运行UMAP需要先安装'umap-learn'包，这里不做介绍，两种方法都可以使用，但不要混用，如果混用，后面的结算结果会将先前的聚类覆盖掉，只能保留一个。
YDL@meta.data$orig.ident<-factor(YDL@meta.data$orig.ident,levels = c("pro","baso"))
##这里采用基于TSNE的聚类方法。
YDL <- RunTSNE(object = YDL, dims = 1:pcSelect,check_duplicates = FALSE)                      #TSNE聚类
pdf(file="TSNE_baso.pdf",width=6.5,height=6)
TSNEPlot(object = YDL, pt.size = 2, label = TRUE)    #TSNE可视化

#另一个可视化的方法
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
dev.off()
write.table(YDL$seurat_clusters,file="tsneCluster_baso.txt",quote=F,sep="\t",col.names=F)


##这里采用基于umap的聚类方法
#这里采用基于图论的聚类方法
YDL <- RunUMAP(object = YDL, dims = 1:pcSelect)                      #umap聚类
pdf(file="umap__baso.pdf",width=6.5,height=6)
UMAPPlot(object = YDL, pt.size = 1.5, label = TRUE)    #umap可视化
#另一个可视化的方法
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident")

DimPlot(object=YDL,label = TRUE,reduction="umap")
##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；
dev.off()

##细胞周期归类

YDL<- CellCycleScoring(object = YDL, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

head(x = YDL@meta.data)
pdf(file="CellCycle_baso.pdf",width=6.5,height=6)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)

DimPlot(YDL,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
dev.off()


pdf("细胞周期打分_pro_baso.pdf")
library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","G2M.Score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = G2M.Score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","S.Score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = S.Score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1

DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 1.5)

DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)

DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",split.by = "orig.ident",pt.size = 1.5)

DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",split.by = "orig.ident",pt.size = 1.5)




library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","G2M.Score","orig.ident")) %>% arrange(G2M.Score)

a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = G2M.Score))+geom_point(size = 1)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~orig.ident )

Idents(YDL)<-YDL@meta.data$Phase

DimPlot(YDL,reduction = "tsne",label = TRUE,split.by = "orig.ident",pt.size = 1.5)

G2M<-subset(YDL,idents = c("G2M"))
DimPlot(G2M,reduction = "tsne",label = TRUE,cols = "#0CB702",split.by = "orig.ident",group.by="Phase",pt.size = 1.5)




library(ggplot2)
mydata<- FetchData(G2M,vars = c("tSNE_1","tSNE_2","G2M.Score","orig.ident")) %>% arrange(G2M.Score)

a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = G2M.Score))+geom_point(size = 1)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~orig.ident )



library(ggplot2)
mydata<- FetchData(G2M,vars = c("tSNE_1","tSNE_2","G2M.Score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = G2M.Score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1






dev.off()

Idents(G2M)<-G2M@meta.data$orig.ident
DimPlot(G2M,reduction = "tsne",label = TRUE,split.by = "orig.ident",group.by="Phase",pt.size = 1.5)

DoHeatmap(object = G2M, features = g2m_genes) + NoLegend()


YDL.AVERAGE<-AverageExpression(object = G2M,return.seurat=F)
YDL.AVERAGE<-as.data.frame(YDL.AVERAGE)

write.csv(YDL.AVERAGE,file = "AVERAGE_pro_baso.csv")
write.csv(YDL.AVERAGE,file = "AVERAGE_neutrophil_baso.csv")

data<-YDL.AVERAGE

merge_tf<-c(g2m_genes)
#挑选部分感兴趣
my.regulons <- merge_tf
#删掉所有列上都重复的
newdata<-data[c(my.regulons),]
newdata<-na.omit(newdata)
newdata<-newdata[which(rowSums(newdata) > 0),]
#colnames(newdata)<-c("pro","baso")
colnames(newdata)<-c("neutrophil","baso")
#低值为蓝色，高值为红色，中间值为白色：
#pdf("fig1g.pdf")
pheatmap(newdata,fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))

pheatmap(newdata,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))




# 加载包
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))




pheatmap(newdata,scale = "row",fontsize = 6,filename = "new2.pdf",width = 10,height = 100,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()




g2m<-read.table("D:/analysis/forpublication/baso_merge_new/g2m_pro_baso.txt")
DoHeatmap(object =G2M, features = g2m$V1) + NoLegend()

g2m<-read.table("D:/analysis/forpublication/baso_merge_new/g2m_neutrophil_baso.txt")
DoHeatmap(object =G2M, features = g2m$V1) + NoLegend()




YDL<-readRDS("D:/analysis/forpublication/baso_merge_new/baso_merge_SINGLET.RDS")

library(Hmisc)
g2m_genes<-capitalize(tolower(cc.genes.updated.2019$g2m.genes))


my_levels <- c(3,2,1,0,4,5)

factor(Idents(YDL), levels= my_levels) 
Idents(YDL) <- factor(Idents(YDL), levels= my_levels)
DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5,cols = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3"))
DoHeatmap(object = YDL, features = g2m_genes,group.colors = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) + NoLegend()


g2m<-read.table("D:/analysis/forpublication/baso_merge_new/g2m_baso.txt")
DoHeatmap(object = YDL, features = g2m$V1,group.colors = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) + NoLegend()


YDL.AVERAGE<-AverageExpression(object = YDL,return.seurat=F)
YDL.AVERAGE<-as.data.frame(YDL.AVERAGE)

write.csv(YDL.AVERAGE,file = "AVERAGE_baso.csv")
data<-YDL.AVERAGE

merge_tf<-c(g2m_genes)
#挑选部分感兴趣
my.regulons <- merge_tf
#删掉所有列上都重复的
newdata<-data[c(my.regulons),]
newdata<-na.omit(newdata)
newdata<-newdata[which(rowSums(newdata) > 0),]
colnames(newdata)<-c("baso_cluster3","baso_cluster2","baso_cluster1","baso_cluster0","baso_cluster4","baso_cluster5")
#低值为蓝色，高值为红色，中间值为白色：
#pdf("fig1g.pdf")
pheatmap(newdata,fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))

pheatmap(newdata,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))




# 加载包
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))




pheatmap(newdata,scale = "row",fontsize = 6,filename = "new2.pdf",width = 10,height = 100,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()







dir.create("./neutrophil")
setwd("./neutrophil")

YDL<-readRDS("D:/analysis/forpublication/baso_merge_new/baso_merge_SINGLET.RDS")
wt_ctl_bm1<-readRDS("D:/analysis/erythropoiesis/wt_ctl_bm1.rds")
wt_ctl_bm1@meta.data$orig.ident<-c("neutrophil")
YDL@meta.data$orig.ident<-c("baso")

head(YDL@meta.data)

YDL<-merge(YDL,wt_ctl_bm1)
Idents(YDL) <- YDL@meta.data$orig.ident
##计算每个细胞的线粒体基因转录本数的百分比（%）,使用[[ ]] 操作符存放到metadata中，mit-开头的为线粒体基因
YDL[["percent.mt"]] <- PercentageFeatureSet(object = YDL, pattern = "^mt-")
###展示基因及线粒体百分比（这里将其进行标记并统计其分布频率，"nFeature_RNA"为基因数，"nCount_RNA"为UMI数，"percent.mt"为线粒体占比）
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


head(YDL@meta.data)
summary(YDL@meta.data$nCount_RNA)
summary(YDL@meta.data$nFeature_RNA)
summary(YDL@meta.data$percent.mt)


##过滤细胞：根据上面小提琴图中基因数"nFeature_RNA"和线粒体数"percent.mt"，分别设置过滤参数，这里基因数 200-4000，线粒体百分比为小于 5%，保留gene数大于200小于2500的细胞；目的是去掉空GEMs和1个GEMs包含2个以上细胞的数据；而保留线粒体基因的转录本数低于5%的细胞,为了过滤掉死细胞等低质量的细胞数据。
YDL <- subset(x = YDL, subset = nFeature_RNA > 200 &nFeature_RNA < 4000 &nCount_RNA<40000 &nCount_RNA>1000 & percent.mt < 15)    #对数据进行过滤
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(Idents(YDL))
summary(YDL@meta.data$nCount_RNA)
summary(YDL@meta.data$nFeature_RNA)
summary(YDL@meta.data$percent.mt)
# > summary(YDL@meta.data$nCount_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1062   17655   20865   21095   24848   39922 
# > summary(YDL@meta.data$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 201    1575    1962    1867    2298    3985 
# > summary(YDL@meta.data$percent.mt)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.03366  0.37025  0.58095  0.95798  0.99198 13.64552 


# YDL<-readRDS("neutrophil_baso_SINGLET.RDS")
# saveRDS(YDL,"neutrophil_baso_SINGLET.RDS")

#对数据进行标准化
##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000
YDL <- NormalizeData(object = YDL, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
YDL <- FindVariableFeatures(object = YDL, selection.method = "vst", nfeatures = 2000)
YDL <- ScaleData(object = YDL, features = rownames(YDL))
#线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
YDL=RunPCA(object= YDL,npcs = 20,pc.genes=VariableFeatures(object = YDL))     #PCA分析
ElbowPlot(YDL)#选择top20个PC
pcSelect=20
YDL <- FindNeighbors(object = YDL, dims = 1:pcSelect)                #计算邻接距离
##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
YDL <- FindClusters(object = YDL, resolution = 0.3)                  #对细胞分组,优化标准模块化
##使用Idents（）函数可查看不同细胞的分群；
head(Idents(YDL), 5)

##Seurat提供了几种非线性降维的方法进行数据可视化（在低维空间把相似的细胞聚在一起）,比如UMAP和t-SNE,运行UMAP需要先安装'umap-learn'包，这里不做介绍，两种方法都可以使用，但不要混用，如果混用，后面的结算结果会将先前的聚类覆盖掉，只能保留一个。

##这里采用基于TSNE的聚类方法。
YDL <- RunTSNE(object = YDL, dims = 1:pcSelect,check_duplicates = FALSE)                      #TSNE聚类
pdf(file="TSNE_baso.pdf",width=6.5,height=6)
TSNEPlot(object = YDL, pt.size = 2, label = TRUE)    #TSNE可视化

#另一个可视化的方法
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
dev.off()
write.table(YDL$seurat_clusters,file="tsneCluster_baso.txt",quote=F,sep="\t",col.names=F)


##这里采用基于umap的聚类方法
#这里采用基于图论的聚类方法
YDL <- RunUMAP(object = YDL, dims = 1:pcSelect)                      #umap聚类
pdf(file="umap__baso.pdf",width=6.5,height=6)
UMAPPlot(object = YDL, pt.size = 1.5, label = TRUE)    #umap可视化
#另一个可视化的方法
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident")

DimPlot(object=YDL,label = TRUE,reduction="umap")
##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；
dev.off()

##细胞周期归类

YDL<- CellCycleScoring(object = YDL, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

head(x = YDL@meta.data)
pdf(file="CellCycle_baso.pdf",width=6.5,height=6)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)

DimPlot(YDL,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
dev.off()

G2M@meta.data$orig.ident<-factor(G2M@meta.data$orig.ident,levels = c("neutrophil","baso"))
