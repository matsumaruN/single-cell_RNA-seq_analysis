###DATASET###-------------------------------------------------------------------
setwd("~/Desktop/R_Seurat/AD.scRNA_Seq/")
library(tidyverse)
library(Seurat)
library(openxlsx)
library(hdf5r)
library(ggplot2)
library(HGNChelper)
library(SeuratData)
set.seed(1234)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); 
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain" 

gs_list = gene_sets_prepare(db_, tissue)
scRNAseqData = readRDS(gzcon(url('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/exampleData.RDS'))); 
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

raw_WT <- Read10X_h5("~/Desktop/Downloaded.Data/AD:WT/GSM7021714_C_6m_filtered_feature_bc_matrix.h5")
raw_AD <- Read10X_h5("~/Desktop/Downloaded.Data/AD:WT/GSM7021717_A_6m_filtered_feature_bc_matrix.h5")

raw_WT <- CreateSeuratObject(counts = raw_WT, project = "WT", min.cells = 3, min.features = 200)
raw_WT$type <- "WT"

raw_AD <- CreateSeuratObject(counts = raw_AD, project = "AD", min.cells = 3, min.features = 200)
raw_AD$type <- "AD"

###QUALITY CHECK###-------------------------------------------------------------
raw_WT[["percent.mt"]] <- PercentageFeatureSet(raw_WT, pattern = "^MT-")
raw_AD[["percent.mt"]] <- PercentageFeatureSet(raw_AD, pattern = "^MT-")

VlnPlot(raw_WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(raw_WT, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(raw_WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(raw_AD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(raw_AD, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(raw_AD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

raw_WT <- subset(raw_WT, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)
raw_AD <- subset(raw_AD, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 25)

###NORMALIZATION###-------------------------------------------------------------
raw_WT <- NormalizeData(raw_WT, normalization.method = "LogNormalize", scale.factor = 10000)
raw_AD <- NormalizeData(raw_AD, normalization.method = "LogNormalize", scale.factor = 10000)

raw_WT <- FindVariableFeatures(raw_WT, selection.method = "vst", nfeatures = 3000)
raw_AD <- FindVariableFeatures(raw_AD, selection.method = "vst", nfeatures = 3000)

top10 <- head(VariableFeatures(raw_WT), 10)
plot1 <- VariableFeaturePlot(raw_WT)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

top10 <- head(VariableFeatures(raw_AD), 10)
plot1 <- VariableFeaturePlot(raw_AD)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

DoHeatmap(merge.data, dims = 1:3, cells = 500, balanced = T)
###統合###----------------------------------------------------------------------
merge.anchors <- FindIntegrationAnchors(object.list = list(raw_WT, raw_AD), dims = 1:20)
merge.data <- IntegrateData(anchorset = merge.anchors, dims = 1:20)
DefaultAssay(merge.data) <- "integrated"

###PCA###-----------------------------------------------------------------------
merge.data <- ScaleData(merge.data, verbose = FALSE)
merge.data <- RunPCA(merge.data, npcs = 30, verbose = FALSE)

###CLUTERING###-----------------------------------------------------------------
merge.data <- FindNeighbors(merge.data, reduction = "pca", dims = 1:20)
merge.data <- FindClusters(merge.data, resolution = 0.5)
merge.data <- RunUMAP(merge.data, reduction = "pca", dims = 1:20)

saveRDS(merge.data,"~/Desktop/Save/WT:AD_Save/mergeData_donePCA.rds")
DimPlot(merge.data, reduction = "umap", split.by = "type", label = TRUE)

###ANNOTATION###----------------------------------------------------------------
merge.data <- readRDS("~/Desktop/Save/WT:AD_Save/mergeData_donePCA.rds")
#陽性遺伝子と陰性遺伝子のスコアリング
es.max = sctype_score(scRNAseqData = merge.data[["integrated"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
#meta.dataのseurat_clustersんも情報とes.maxの情報を統合
cL_resutls = do.call("rbind", lapply(unique(merge.data@meta.data$seurat_clusters), function(cl){
  es.max.cl = 
sort(rowSums(es.max[ ,rownames(merge.data@meta.data[merge.data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(merge.data@meta.data$seurat_clusters==cl)), 10)
}))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 
#信頼性の低いsctype_scoreの場合そのクラスターをunknownに
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < 
sctype_scores$ncells/4] = "unknown"
  print(sctype_scores[,1:3])
#細胞腫を特定する新たなmeta.dataである"customclassif"を追加する
merge.data@meta.data$customclassif = ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j,]
  merge.data@meta.data$customclassif[merge.data@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}

DimPlot(merge.data, reduction = "umap", split.by = "type", label = TRUE, repel = TRUE, group.by = 'customclassif') 

###アノテーションした細胞腫を抜き出して、DEGを見つける###-----------------------
Idents(merge.data) <- "customclassif"
#merge.data2 <- subset(merge.data, idents = 0)
#Dimplotで確認できる
MC_data <- subset(merge.data, idents = "Microglial cells")
#MC_data <- SetIdent(MC_data, value = "type")
MC_list<- SplitObject(object = MC_data, split.by = "type")
Galectin <- subset(merge.data, idents = "LGALS3")
#DEGの探索とp値のFDRによるBH法補正
MC_marker <- FindMarkers(MC_data, ident.1 = "WT", assay = "RNA", lfcThreshold = 0, min.cells.feature = 1, min.cells.group = 1, min.pct = 0)
MC_marker$p_val_adj = p.adjust(MC_marker$p_val, method = 'fdr')
write.csv(MC_marker, "RNAMC_DEG.csv")
VlnPlot(MC_data, "Galectin3", assay = "RNA")

DoHeatmap(object = MC_data, group.by = "type", slot = "scale.data")
###merge.dataをWTとADで分割してlist化する###------------------------------------
type_data <- SplitObject(merge.data, split.by = "type")
#tmpとはtemporaryのこと、意味のない一時的なファイル
#table()を使うと各インデントの回数を"Freq"として保存。
#Var1,Var2列も自動的で生成される。それぞれ"ident"と"type"が格納されている
tmp1 <- as.data.frame(table(Idents(type_data$WT), type_data$WT$type))
tmp2 <- as.data.frame(table(Idents(type_data$AD), type_data$AD$type))
#tmpをbind_rows()で結合させる
pt <- bind_rows(tmp1, tmp2)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank())

MC_data <- subset(merge.data, idents = "Microglial cells")
MC_data <- SetIdent(MC_data, value = "type")

###各クラスターのDEGのエンリッチ解析と比較###-----------------------------------
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

deg <- FindAllMarkers(object = merge.data, only.pos = T)
table(deg$cluster)
Seurat::DoHeatmap(raw_AD, features = gene_list,)
#クラスターごとのDEGをベクトルにしてリストにまとめる
gene_list <- split(deg$gene, deg$cluster)
#遺伝子名から遺伝子IDに変換
for(i in names(gene_list)){
  gene_list[[i]] <- bitr(gene_list[[i]],
                           fromType = "SYMBOL",
                           toType = "ENTREZID",
                           OrgDb = "org.Hs.eg.db")
                           }
#残った遺伝子を調べる
length(unlist(gene_list))/2

#clusterProfiler用に成形
for(i in names(gene_list)){
  gene_list[[i]] <- gene_list[[i]]$ENTREZID
}

#Gene Ontology
GOBP <- compareCluster(geneClusters = gene_list,
                          fun = "enrichGO",
                          ont = "BP",
                          OrgDb = "org.Hs.eg.db")

#KEGG
kegg <- compareCluster(geneClusters = gene_list,
                       fun = "enrichKEGG",
                       organism = "mmu")

#Reactome
library(ReactomePA)
reactome <- compareCluster(geneClusters = gene_list,
                           fun = "enrichPathway",
                           organism = "mouse")

#WikiPathway
WP_clusterplot <- compareCluster(geneClusters = gene_list,
                                 fun = "enrichWP",
                                 organism = "Mus musculus")

#Disease Ontology
library(DOSE)
DO <- compareCluster(geneClusters = gene_list,
                     fun = "enrichDO")
#DisGeNET
DGN <- compareCluster(geneClusters = gene_list,
                      fun = "enrichDGN")
#結果
dotplot(GOBP)


###volcanoplot###---------------------------------------------------------------
result <- readRDS("~/Desktop/Save/WT:AD_Save/mergeData_donePCA.rds")
selectLab <- result %>% filter(PValue < 0.05) %>% filter(abs(logFC) > 2.5) %>% .$extemal_gene_name


FeaturePlot(object = merge.data, features = "Acta2", label = TRUE, reduction = "umap", split.by = "type", group.by = 'customclassif')

