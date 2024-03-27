library(scran)
library(irlba)
library(ggplot2)
library(batchelor)

library(umap)
library(Seurat)
library(reticulate)
use_condaenv(condaenv="scanpy-p3.9")

umap = import('umap')

###
path2data1  <- "/data1/ivanir/Chp2022/10xGen/analysis/"

sce_10x <- readRDS(paste0(path2data1, "sce.rds"))
sce_10x <- sce_10x[, colData(sce_10x)$scDblFinder.class == "singlet"]

path2data2  <- '/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/Combined/all-well/DGE_filtered/'

sce_parse <- readRDS(paste0(path2data2, "sce.rds"))
sce_parse <- sce_parse[, colData(sce_parse)$scDblFinder.class == "singlet"]

genes_shared <- intersect(rownames(sce_10x), rowData(sce_parse)$gene_name)
sce_10x   <- sce_10x[genes_shared, ]
sce_parse <- sce_parse[match(genes_shared,rowData(sce_parse)$gene_name), ]
rownames(sce_parse) <- genes_shared
 
#prevent duplicate rownames
colnames(sce_10x)   <- paste0("10X_", colnames(sce_10x))
colnames(sce_parse) <- paste0("PB_",  colnames(sce_parse))
  
sce <- SingleCellExperiment::SingleCellExperiment(
      list(counts=cbind(counts(sce_10x), counts(sce_parse))))

meta_10x   <- data.frame(colData(sce_10x))[,c("sequencing.round","sample","day","mt.fraction","scDblFinder.score")]
meta_parse <- data.frame(colData(sce_parse))[,c("seq_library","sample_name","day","mt.fraction","scDblFinder.score")]

colnames(meta_10x) <- colnames(meta_parse) <- c("batch","sample","day","mt.fraction","doublet.score")

meta <- rbind(meta_10x,meta_parse)
lib.sizes <- colSums(counts(sce))
ngenes    <- colSums(counts(sce) > 0)

platform  <- rep("10X", nrow(meta))
platform[grep("DNA", meta$batch)] <- "PB"
meta <- cbind(meta, ngenes, lib.sizes, platform)

rownames(meta) <- colnames(sce)
colData(sce)   <- DataFrame(meta)

rowData(sce) <- rowData(sce_parse)

saveRDS(sce,"/data1/ivanir/Chp2022/Integrated/sce.rds")

writeLines(colnames(sce),"/data1/ivanir/Chp2022/Integrated/cell_names.txt")
writeLines(rownames(sce),"/data1/ivanir/Chp2022/Integrated/gene_names.txt")

write.table(data.frame(colData(sce)), "/data1/ivanir/Chp2022/Integrated/meta.tab",
		quote=FALSE, row.names=FALSE, sep="\t")

sce <- logNormCounts(sce)

write.table(as.matrix(logcounts(sce)), "/data1/ivanir/Chp2022/Integrated/normalised_counts.tab",
		quote=FALSE, row.names=FALSE, sep="\t")

write.table(as.matrix(counts(sce)), "/data1/ivanir/Chp2022/Integrated/counts.tab",
		quote=FALSE, row.names=FALSE, sep="\t")

###
###
sce  <- readRDS("/data1/ivanir/Chp2022/Integrated/sce.rds")
meta <- data.frame(colData(sce))

pdf("/data1/ivanir/Chp2022/Integrated/UMIsByPlatform.pdf")
ggplot(meta, aes (x = factor(platform), y = as.numeric(lib.sizes))) +
  geom_boxplot() +
  theme_bw() + 
  labs(x = "Platform", y = "Number of UMIs") +
  scale_y_log10(breaks = c(1000, 5000, 10000, 50000, 100000),
    labels = c("1,000", "5,000", "10,000", "50,000", "100,000"))
dev.off()

pdf("/data1/ivanir/Chp2022/Integrated/MitFraction.pdf")
ggplot(meta, aes (x = factor(platform), y = as.numeric(mt.fraction))) +
  geom_boxplot() +
  theme_bw() + 
  labs(x = "Platform", y = "Mitochondrial fraction") 
dev.off()

pdf("/data1/ivanir/Chp2022/Integrated/cell_complexity.pdf")
Platform <- factor(meta$platform)
qplot(meta$lib.sizes, meta$ngenes, col=Platform) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() + 
  theme(text = element_text(size=20))  +
  labs(x = "UMI count", y = "Number of expressed genes")
dev.off()

sce <- sce[calculateAverage(sce)>0.05,]
sce <- logNormCounts(sce)

decomp  <- modelGeneVar(sce)
hvgs    <- rownames(decomp)[decomp$FDR < 0.05]
pca     <- prcomp_irlba(t(logcounts(sce[hvgs,])), n = 30)
layout  <- umap(pca$x, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs", "min_dist"), n_neighbors=15, min_dist=.99)

df_plot <- data.frame(
 meta,
 UMAP1 = layout$layout[,1],
 UMAP2 = layout$layout[,2] 
)

plot.index <- sample(nrow(df_plot))

ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = batch)) +
  geom_point(size = 0.4) +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("/data1/ivanir/Chp2022/Integrated/umap_by_batch.pdf")

day_colours <- rev(wesanderson::wes_palette("Zissou1", 15, type = "continuous"))
days <- sort(as.numeric(unique(meta$day)))

level_order <- days
df_plot$day <- factor(df_plot$day,levels=level_order)

ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = day)) +
  geom_point(size = 0.4) +
  labs(x = "Dim 1", y = "Dim 2") +
  scale_color_manual(values=day_colours, name = "Day") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("/data1/ivanir/Chp2022/Integrated/umap_by_day.pdf")

df_plot <- data.frame(
  day   = names(table(meta$day)),
  n_cells = as.numeric(as.character(table(meta$day)))
)

level_order <- days
df_plot$day <- factor(df_plot$day,levels=level_order)

ggplot(df_plot, aes(x=day , y=n_cells, fill=day)) +
  geom_bar(stat="identity")+theme_minimal() +
  scale_fill_manual(values = day_colours,name = "Day") +
  labs(x = "Day", y = "Number of cells")
ggsave("/data1/ivanir/Chp2022/Integrated/cells_by_day.pdf")

## MNN batch correction
cos_normalised <- sapply(sort(as.numeric(as.character(unique(df_plot$day)))), 
  function(x) {sce_sub <- logcounts(sce[hvgs,which(df_plot$day == x)]); cosineNorm(sce_sub)})
names(cos_normalised)  <- sort(as.numeric(as.character(unique(df_plot$day))))
pcs <- multiBatchPCA(cos_normalised)

mnn.out <- reducedMNN(pcs, merge.order = 15:1)

layout  <- umap(mnn.out$corrected, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs", "min_dist"), n_neighbors=15, min_dist=.99)

df_plot <- data.frame(
 meta,
 UMAP1 = layout$layout[,1],
 UMAP2 = layout$layout[,2] 
)

ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = batch)) +
  geom_point(size = 0.4) +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("/data1/ivanir/Chp2022/Integrated/umap_by_batch_mnn.pdf")

day_colours <- rev(wesanderson::wes_palette("Zissou1", 15, type = "continuous"))
days <- sort(unique(meta$day))

level_order <- days
df_plot$day <- factor(df_plot$day,levels=level_order)

ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = day)) +
  geom_point(size = 0.4) +
  labs(x = "Dim 1", y = "Dim 2") +
  scale_color_manual(values=day_colours, name = "Day") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("/data1/ivanir/Chp2022/Integrated/umap_by_day_mnn.pdf")

#mnn.out.exp <- mnnCorrect(sce, batch=meta$batch)
#dim(mnn.out.exp@assays@data$corrected)

### Seurat batch correction
seurat_integ <- CreateSeuratObject(as.matrix(counts(sce)), meta.data = data.frame(meta))
seurat_list_platform <- SplitObject(seurat_integ, split.by = "platform")[2]
seurat_list_sample   <- SplitObject(seurat_integ, split.by = "sample")[1:9]
seurat_list          <- c(seurat_list_sample,seurat_list_platform)
#seurat_list <- scCustomize::Merge_Seurat_List(list_seurat = c(seurat_list_sample,seurat_list_platform), merge = FALSE)

# normalize and identify variable features for each dataset independently
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = seurat_list)
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

seurat_list_10X <- seurat_list[1:9]
seurat_list_PB  <- seurat_list[10]
 
anchors <- FindIntegrationAnchors(object.list = seurat_list_10X, anchor.features = features, reduction = "rpca")

combined_10x <- IntegrateData(anchorset = anchors)
DefaultAssay(combined_10x) <- "integrated"

seurat_list <- c(combined_10x,seurat_list_PB)
names(seurat_list) <- c("10X","PB")
features <- SelectIntegrationFeatures(object.list = seurat_list)
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

#anchors  <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, reduction = "rpca")
anchors  <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, reduction = "cca")

combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

pca_corrected <- combined@reductions$pca@cell.embeddings
colnames(pca_corrected) <- paste0("PC",1:30)
write.table(pca_corrected, file="/data1/ivanir/Chp2022/Integrated/pca_corrected_seurat_cca.csv", row.names=FALSE, col.names=TRUE, sep=",")
saveRDS(combined,"/data1/ivanir/Chp2022/Integrated/seurat_object_cca.rds")

layout  <- umap(pca_corrected, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs", "min_dist"), n_neighbors=15, min_dist=.99)

df_plot <- data.frame(
 meta,
 UMAP1 = layout$layout[,1],
 UMAP2 = layout$layout[,2] 
)
plot.index <- sample(nrow(df_plot))

ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = batch)) +
  geom_point(size = 0.4) +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("/data1/ivanir/Chp2022/Integrated/umap_by_batch_seurat_cca.pdf")

day_colours <- rev(wesanderson::wes_palette("Zissou1", 15, type = "continuous"))
days <- sort(as.numeric(unique(meta$day)))

level_order <- days
df_plot$day <- factor(df_plot$day,levels=level_order)

ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = day)) +
  geom_point(size = 0.4) +
  labs(x = "Dim 1", y = "Dim 2") +
  scale_color_manual(values=day_colours, name = "Day") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("/data1/ivanir/Chp2022/Integrated/umap_by_day_seurat_cca.pdf")

### Regress out cell cycle
library(Seurat)
s_genes   <- c('MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8')
g2m_genes <- c('HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA')

combined <- readRDS("/data1/ivanir/Chp2022/Integrated/seurat_object_cca.rds")

combined <- CellCycleScoring(combined, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)

combined <- ScaleData(combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(combined))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
combined <- RunPCA(combined, features = VariableFeatures(combined), npcs = 30, nfeatures.print = 10)

pca_corrected <- combined@reductions$pca@cell.embeddings
colnames(pca_corrected) <- paste0("PC",1:30)
write.table(pca_corrected, file="/data1/ivanir/Chp2022/Integrated/pca_corrected_seurat_cca_cellcycle_reg.csv", row.names=FALSE, col.names=TRUE, sep=",")
saveRDS(combined,"/data1/ivanir/Chp2022/Integrated/seurat_object_cca_cellcycle_reg.rds")
