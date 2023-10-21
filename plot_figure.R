
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

TSNEPlot(object = YDL, pt.size = 2, label = TRUE)    #TSNE可视化



library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
DimPlot(YDL, reduction = "tsne",pt.size = 1.5,label = T)

DimPlot(YDL, reduction = "umap",pt.size = 1.5,label = T)
###Pseudotime monocle3
cds <- as.cell_data_set(YDL)
cds <- cluster_cells(cds)
head(pData(cds))
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)


stem1<-rownames(pData(cds)[which(pData(cds)$ident %in% c('3')),])
cds <- order_cells(cds, root_cells = stem1)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
pdf("monocle3-baso.pdf")
dev.off()



YDL<-subset(YDL,idents = c("0","1","2","3","4","5"))
a<-as.data.frame(table(YDL@meta.data$seurat_clusters))
a<-a[-7,]
# 数据准备
info = a$Freq
# 命名
#names = as.character(a$Var1)
names = c("0","1","2","3","4","5")
# 涂色（可选）
#cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC","#3f90CC","#22B17C","#FFC93E")
# 计算百分比
piepercent = paste(round(100*info/sum(info),2), "%")

# 绘图
pie(info, labels=piepercent, main = "total cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.6, fill=colors)
pdf("figure5A.pdf")
DimPlot(YDL, reduction = "tsne",pt.size = 1.5,label = T)
# 绘图
pie(info, labels=piepercent, main = "total cluster ratio", col=colors, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.6, fill=colors)

dev.off()

pdf("figure5B.pdf",height = 12,width = 24)

FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Hbb-bt","Hbb-bs","Bpgm","Alas2","Tcp11l2","Cd36","Xpo7","Tbcel"),cols = c("gray", "red"),ncol = 4)#actin

dev.off()



pdf("figure7B.pdf",height = 12,width = 18)

FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Hist1h1b","Hist1h2ae","Hist1h2bn","Hist1h3c","Hist1h4d","Fbxo5"),cols = c("gray", "red"),ncol = 3)#actin

dev.off()



library(monocle)
cds<-readRDS("cds_SINGLET.rds")

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


pdf("figure5C.pdf")
plot_cell_trajectory(cds, color_by = "Cluster")
plot_cell_trajectory(cds, color_by="Pseudotime", show_backbone=FALSE)
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




pdf("figure5D.pdf")
DoHeatmap(object = subset(YDL, downsample = 100), features = top10$gene,group.colors = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) + NoLegend()
dev.off()

gene<- YDL.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

DotPlot(object = YDL, features = unique(gene$gene))+ RotatedAxis()
DotPlot(object = YDL, features =unique(gene$gene))+ RotatedAxis()
DotPlot(object = YDL, features =unique(gene$gene))+ RotatedAxis()+ coord_flip()

#加边框
DotPlot(object = YDL, features =unique(gene$gene))+ coord_flip()+ RotatedAxis()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)


pdf("figure5E.pdf")
DotPlot(object = YDL, features =unique(gene$gene))+ coord_flip()+ RotatedAxis()+annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+ annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 0.5)
dev.off()




pdf("figure5F.pdf")

DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)
#DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")

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
dev.off()



pdf("figure5H.pdf")
library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","G2M.Score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = G2M.Score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1
dev.off()







library(Hmisc)
g2m_genes<-capitalize(tolower(cc.genes.updated.2019$g2m.genes))

pdf("figure5I.pdf")
DoHeatmap(object = YDL, features = g2m_genes,group.colors = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) + NoLegend()
DoHeatmap(object = subset(YDL, downsample = 100), features = g2m_genes,group.colors = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")) + NoLegend()
dev.off()




YDL<-readRDS("D:/analysis/forpublication/baso_merge_new/pro/neutrophil/neutrophil_baso_SINGLET.RDS")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)

DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE,group.by = 'orig.ident')
DimPlot(YDL, reduction = "tsne", pt.size = 1.5, label=TRUE,group.by = 'orig.ident')




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

ALL<-readRDS("D:/analysis/forpublication/bulk_sc/terminalE.rds")
Idents(ALL)<-ALL$orig.ident
pro<-subset(ALL,idents = c("pro"))
pro@meta.data$orig.ident<-c("pro")
rm(ALL)
wt_ctl_bm1<-readRDS("D:/analysis/erythropoiesis/wt_ctl_bm1.rds")
wt_ctl_bm1@meta.data$orig.ident<-c("neutrophil")
YDL@meta.data$orig.ident<-c("baso")

head(YDL@meta.data)

YDL<-merge(YDL,wt_ctl_bm1)
YDL<-merge(YDL,pro)
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
YDL <- subset(x = YDL, subset = nFeature_RNA > 200 &nFeature_RNA < 8000 &nCount_RNA<100000 &nCount_RNA>1000 & percent.mt < 15)    #对数据进行过滤
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


# YDL<-readRDS("PRO_BASO_NEUTROPHIL.RDS")
# saveRDS(YDL,"baso_SINGLET.RDS")

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

DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "orig.ident")

DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "Phase")

pdf("figure5G-1.pdf")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident",group.by = "Phase")
dev.off()


Idents(YDL)<-YDL@meta.data$Phase

DimPlot(YDL,reduction = "tsne",label = TRUE,split.by = "orig.ident",pt.size = 1.5)

G2M<-subset(YDL,idents = c("G2M"))



pdf("figure5G-2.pdf")
DimPlot(G2M,reduction = "tsne",label = TRUE,cols = "#0CB702",split.by = "orig.ident",group.by="Phase",pt.size = 1.5)
dev.off()


YDL@meta.data$sample<-YDL@meta.data$orig.ident
library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","G2M.Score","orig.ident")) %>% arrange(G2M.Score)

a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = G2M.Score))+geom_point(size = 1)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

pdf("figure5G-3.pdf")
a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~orig.ident) +
  theme(
    strip.background = element_rect(
      color = "white", fill = "white"),
    panel.grid = element_blank())
dev.off()





library(Hmisc)
g2m_genes<-capitalize(tolower(cc.genes.updated.2019$g2m.genes))


setwd("D:/analysis/forpublication/baso_merge_new")
YDL<-readRDS("baso_merge_SINGLET.RDS")
my_levels <- c(3,2,1,0,4,5)
factor(Idents(YDL), levels= my_levels) 
Idents(YDL) <- factor(Idents(YDL), levels= my_levels)
DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5,cols = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3"))

YDL.AVERAGE<-AverageExpression(object = YDL,return.seurat=F)
YDL.AVERAGE<-as.data.frame(YDL.AVERAGE)



write.csv(YDL.AVERAGE,file = "AVERAGE_baso.csv")
data<-YDL.AVERAGE

data<-data[which(rowSums(data) > 0),]#去掉全为零的行 情况

merge_tf<-c(g2m_genes)
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

pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))




pheatmap(newdata,scale = "row",fontsize = 6,filename = "new2.pdf",width = 10,height = 100,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()

pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))


type1 <- c(
  "Psrc1", "Hjurp", "Cbx5", "Ckap2l", "Kif11", "Ndc80", "Rangap1", "Smc4", 
  "Bub1", "Birc5", "Mki67", "Hmgb2", "Anp32e", "Cdca3", "Aurkb", "Ttk", 
  "Top2a", "Ncapd2", "Ctcf", "Tmpo", "Lbr"
)

type2 <- c(
  "Ccnb2", "Dlgap5", "Tubb4b", "Nek2", "Jpt1", "Gas2l3", "Ect2", "Cenpa", 
  "G2e3", "Ube2c", "Cenpf", "Cks1b", "Cenpe", "Cdc20", "Nusap1", "Kif2c", 
  "Cdk1", "Gtse1", "Cdca2", "Kif23", "Anln", "Tpx2", "Cks2", "Cdca8", 
  "Ckap5", "Nuf2", "Tacc3", "Cdc25c", "Kif20b", "Ckap2", "Hmmr", "Aurka"
)

eg = bitr(type1, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
type1 <- eg[,2]
type1<-as.data.frame(type1)
colnames(type1)<-c("type1")
head(type1)
type1<-c(type1)


eg = bitr(type2, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
type2 <- eg[,2]
type2<-as.data.frame(type2)
colnames(type2)<-c("type2")
head(type2)
type2<-c(type2)


data<-list(type1 =type1$type1,type2=type2$type2)

lapply(data, head)


head(as.data.frame(ck))

ck <- compareCluster(geneCluster = data,OrgDb = org.Mm.eg.db, fun = "enrichGO", pvalueCutoff=0.05)
write.csv(ck,"go.csv")


dim(ck)
pdf("figure5K.pdf",height = 7,width=8)
dotplot(ck, showCategory =50)
dev.off()
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





























































YDL<-readRDS("D:/analysis/forpublication/bulk_sc/terminalE.rds")
#徐湘民的骨髓单细胞数据
#YDL<-readRDS("D:/analysis/forpublication/summary/nonpro1_own_nonpro2/terminalE.rds")
#saveRDS(YDL,"terminalE.rds")
Idents(YDL)<-YDL$orig.ident
YDL<-subset(YDL,idents = c("pro","baso","poly","ortho"))

Idents(YDL)<-YDL$celltype
YDL<-subset(YDL,idents = c("E0","E1","E2","E3","E4"))
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE)
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE,group.by = 'orig.ident')


cols = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE,cols = cols)

pdf("去除reti.pdf")
p1<-DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE)
p1
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE,cols = cols)
dev.off()

pdf("figure4c-1.pdf")
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE,cols = cols)
Idents(YDL)<-YDL$orig.ident
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE,cols = cols,group.by = 'orig.ident')
dev.off()
pdf("figure4c-2.pdf",height = 6,width = 24)
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE,cols = cols,split.by = 'orig.ident')
dev.off()


pdf("figur6-IJK.pdf")
Idents(YDL)<-YDL$orig.ident
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE,cols = cols,group.by = 'orig.ident')
Idents(YDL)<-YDL$celltype
DimPlot(YDL, reduction = "umap", pt.size = 1.5, label=TRUE,cols = cols)
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Atm"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Atr"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Foxm1"),cols = c("gray", "red"))#actin
dev.off()


# 
# library("Seurat")
# library("ggplot2")
# gg <- TSNEPlot(YDL)
# col<-ggplot_build(gg)$data
# col<-as.data.frame(col)
# table(col$colour)
# table(YDL$celltype)
# > table(col$colour)
# 
# #00BA38 #00BFC4 #619CFF #B79F00 #F564E3 #F8766D 
# 2138    2050    3189    1037     696    2564 
# > table(YDL$celltype)
# 
# E0   E1   E2   E3   E4   E5 
# 2564 1037 2138 2050 3189  696 



library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nCount_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nFeature_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = log(nFeature_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = nCount_RNA))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nFeature_RNA"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = nFeature_RNA))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




pdf("Figure_S3_D_E.pdf")


library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nFeature_RNA"))
head(mydata)
colnames(mydata)<-c("UMAP_1","UMAP_2","gene_number")
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = gene_number))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
library(ggplot2)
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2","nCount_RNA"))
head(mydata)
colnames(mydata)<-c("UMAP_1","UMAP_2","mRNA_number")
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = mRNA_number))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

dev.off()
















