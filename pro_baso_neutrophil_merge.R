#




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


setwd("D:/analysis/forpublication/baso_merge_new/merge")


YDL<-readRDS("D:/analysis/forpublication/baso_merge_new/pro/pro_baso_SINGLET.RDS")
wt_ctl_bm1<-readRDS("D:/analysis/erythropoiesis/wt_ctl_bm1.rds")
wt_ctl_bm1@meta.data$orig.ident<-c("neutrophil")
#YDL@meta.data$orig.ident<-c("baso")

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





head(YDL@meta.data)

write.csv(YDL@meta.data,"细胞周期打分-baso+pro+neutrophil-G2M期细胞.csv")
write.csv(YDL@meta.data,"细胞周期打分-baso+pro+neutrophil-所有细胞.csv")



pro<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="pro",]))

baso<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="baso",]))

neutrophil<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="neutrophil",]))


table(YDL@meta.data$Phase)
table(pro@meta.data$Phase)
table(baso@meta.data$Phase)
table(neutrophil@meta.data$Phase)



pdf("细胞周期比例-各个样本.pdf")
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
pie(info, labels=piepercent, main = "all cell cycle ratio", col=cols, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.8, fill=cols)



a<-as.data.frame(table(pro@meta.data$Phase))
# 数据准备
info = a$Freq
# 命名
names = c("G1","G2M","S")
# 涂色（可选）
cols = c("#F8766D","#0CB702","#00A9FF")
# 计算百分比
piepercent = paste(round(100*info/sum(info),2), "%")
# 绘图
pie(info, labels=piepercent, main = "pro cell cycle ratio", col=cols, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.8, fill=cols)




a<-as.data.frame(table(baso@meta.data$Phase))
# 数据准备
info = a$Freq
# 命名
names = c("G1","G2M","S")
# 涂色（可选）
cols = c("#F8766D","#0CB702","#00A9FF")
# 计算百分比
piepercent = paste(round(100*info/sum(info),2), "%")
# 绘图
pie(info, labels=piepercent, main = "baso cell cycle ratio", col=cols, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.8, fill=cols)




a<-as.data.frame(table(neutrophil@meta.data$Phase))
# 数据准备
info = a$Freq
# 命名
names = c("G1","G2M","S")
# 涂色（可选）
cols = c("#F8766D","#0CB702","#00A9FF")
# 计算百分比
piepercent = paste(round(100*info/sum(info),2), "%")
# 绘图
pie(info, labels=piepercent, main = "neutrophil cell cycle ratio", col=cols, family='GB1')
# 添加颜色样本标注
legend("topright", names, cex=0.8, fill=cols)





dev.off()




Idents(YDL)<-YDL@meta.data$Phase
DimPlot(YDL,reduction = "tsne",label = TRUE,split.by = "orig.ident",pt.size = 1.5)

YDL<-subset(YDL,idents = c("G2M"))
DimPlot(YDL,reduction = "tsne",label = TRUE,cols = "#0CB702",split.by = "orig.ident",group.by="Phase",pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,cols = "#0CB702",group.by="Phase",pt.size = 1.5)


table(YDL@meta.data$Phase)
summary(YDL@meta.data$G2M.Score)
summary(YDL@meta.data$S.Score)
summary(YDL@meta.data[YDL@meta.data$Phase=="G2M",]$G2M.Score)
summary(YDL@meta.data[YDL@meta.data$Phase=="S",]$S.Score)
summary(YDL@meta.data[YDL@meta.data$Phase=="G2M",]$S.Score)
summary(YDL@meta.data[YDL@meta.data$Phase=="S",]$G2M.Score)



summary(pro@meta.data$G2M.Score)

summary(baso@meta.data$G2M.Score)

summary(neutrophil@meta.data$G2M.Score)


#Baso的G2M分大于1和大于0.5的占比

table(YDL@meta.data$G2M.Score>1)
304/(304+8070)*100#百分比
table(YDL@meta.data$G2M.Score>0.5)
3338/(3338+5036)*100#百分比

table(YDL@meta.data[YDL@meta.data$Phase=="G2M",]$G2M.Score>1)
304/(304+3799)*100#百分比
table(YDL@meta.data[YDL@meta.data$Phase=="G2M",]$G2M.Score>0.5)
3338/(3338+765)*100#百分比



#Baso的G2M分大于1和大于0.5的占比

table(neutrophil@meta.data$G2M.Score>0.6)
2/(2+175)*100#百分比
table(pro@meta.data$G2M.Score>0.6)
10/(10+93)*100#百分比
table(baso@meta.data$G2M.Score>0.6)
2732/(2732+1144)*100#百分比



#Baso的G2M分大于1和大于0.5的占比

table(YDL@meta.data$G2M.Score>1)
304/(304+8070)*100#百分比
table(YDL@meta.data$G2M.Score>0.5)
3338/(3338+5036)*100#百分比

table(YDL@meta.data[YDL@meta.data$Phase=="G2M",]$G2M.Score>1)
304/(304+3799)*100#百分比
table(YDL@meta.data[YDL@meta.data$Phase=="G2M",]$G2M.Score>0.5)
3338/(3338+765)*100#百分比



pro<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="pro",]))

baso<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="baso",]))


neutrophil<-subset(YDL, cells= rownames(YDL@meta.data[YDL@meta.data$orig.ident=="neutrophil",]))



data<-as.data.frame(YDL@meta.data$G2M.Score)
colnames(data)<-"G2M.score"
library(ggplot2)
#注释：package使用之前需要调用

p<-ggplot(data, aes(G2M.score)) +
  geom_histogram(breaks=seq(-0.7,1.5,0.1))+ xlim(-0.7,1.5)
p


data<-as.data.frame(pro@meta.data$G2M.Score)
colnames(data)<-"G2M.score"
library(ggplot2)
#注释：package使用之前需要调用

p<-ggplot(data, aes(G2M.score)) +
  geom_histogram(breaks=seq(-0.7,1.5,0.1))+ xlim(-0.7,1.5)
p


data<-as.data.frame(baso@meta.data$G2M.Score)
colnames(data)<-"G2M.score"
library(ggplot2)
#注释：package使用之前需要调用

p<-ggplot(data, aes(G2M.score)) +
  geom_histogram(breaks=seq(-0.7,1.5,0.1))+ xlim(-0.7,1.5)
p


data<-as.data.frame(neutrophil@meta.data$G2M.Score)
colnames(data)<-"G2M.score"
library(ggplot2)
#注释：package使用之前需要调用

p<-ggplot(data, aes(G2M.score)) +
  geom_histogram(breaks=seq(-0.7,1.5,0.1))+ xlim(-0.7,1.5)
p


saveRDS(YDL,"PRO_BASO_NEUTROPHIL.RDS")



library(Hmisc)
g2m_genes<-capitalize(tolower(cc.genes.updated.2019$g2m.genes))

DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)

Idents(YDL) <- YDL$orig.ident
DoHeatmap(object = YDL, features = g2m_genes) + NoLegend()


g2m<-read.table("g2m.txt")
DoHeatmap(object = YDL, features = g2m$V1) + NoLegend()

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
colnames(newdata)<-c("baso","pro","neutrophil")
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


pheatmap(newdata,scale = "row",fontsize = 7,#clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))





pdf("细胞周期打分.pdf")
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


Idents(YDL)<-YDL@meta.data$Phase

DimPlot(YDL,reduction = "tsne",label = TRUE,split.by = "orig.ident",pt.size = 1.5)

G2M<-subset(YDL,idents = c("G2M"))
DimPlot(G2M,reduction = "tsne",label = TRUE,cols = "#0CB702",split.by = "orig.ident",group.by="Phase",pt.size = 1.5)

library(ggplot2)
mydata<- FetchData(G2M,vars = c("tSNE_1","tSNE_2","G2M.Score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = G2M.Score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1


baso1<-subset(G2M, cells= rownames(G2M@meta.data[G2M@meta.data$orig.ident=="baso1",]))
baso2<-subset(G2M, cells= rownames(G2M@meta.data[G2M@meta.data$orig.ident=="baso2",]))
neutrophil<-subset(G2M, cells= rownames(G2M@meta.data[G2M@meta.data$orig.ident=="neutrophil",]))


library(ggplot2)
mydata<- FetchData(neutrophil,vars = c("tSNE_1","tSNE_2","G2M.Score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = G2M.Score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1

library(ggplot2)
mydata<- FetchData(baso1,vars = c("tSNE_1","tSNE_2","G2M.Score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = G2M.Score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1


library(ggplot2)
mydata<- FetchData(baso2,vars = c("tSNE_1","tSNE_2","G2M.Score"))
a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = G2M.Score))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

p1<-a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2<-DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
p1




dev.off()


table(YDL@meta.data$orig.ident)
YDL@meta.data[YDL@meta.data$orig.ident=="baso1"]<-"baso"
YDL@meta.data[YDL@meta.data$orig.ident=="baso2"]<-"baso"


Idents(YDL)<-YDL@meta.data$Phase

DimPlot(YDL,reduction = "tsne",label = TRUE,split.by = "orig.ident",pt.size = 1.5)

YDL<-subset(YDL,idents = c("G2M"))
YDL@meta.data$sample<-YDL@meta.data$orig.ident
library(ggplot2)
mydata<- FetchData(YDL,vars = c("tSNE_1","tSNE_2","G2M.Score","orig.ident")) %>% arrange(G2M.Score)

a <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = G2M.Score))+geom_point(size = 1)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~orig.ident )

dev.off()

