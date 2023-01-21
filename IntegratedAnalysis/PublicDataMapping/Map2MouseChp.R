library(scran)
library(irlba)
library(Matrix)
library(SingleCellExperiment)

#library("leiden")
#################################
### Low level analysis
#################################
 
# Read metadata

setwd('/data1/ivanir/Chp2022/ChpMouseCell2021')
counts <-  as.matrix(readMM("sc_rna_counts.mtx"))
cells <- readLines("sc_rna_barcodes.tsv")
genes <- read.table("sc_rna_counts_genes.tsv")
rownames(counts) <- genes$V1
colnames(counts) <- cells

meta_Mchp <- readLines("scEmbryo_cell_metadata_edited.txt")
cells_index <- grep(".Rep", meta_Mchp)
cells       <- meta_Mchp[cells_index]
celltypes   <- meta_Mchp[cells_index+3]

cell_names <- unlist(lapply(strsplit(cells, split="\t"), function(x){x[1]}))

celltypes[grep("epithelial cell",celltypes)] <- "Epithelial cell"
celltypes[grep("endothelial cell",celltypes)] <- "Endothelial cell"                                             
celltypes[grep("mesenchymal",celltypes)] <- "Mesenchymal cell"
celltypes[grep("Oligodendrocyte precursor",celltypes)] <- "Oligodendrocyte precursor"
celltypes[grep("Cycling G2/M",celltypes)] <- "Neuron associated Cycling G2/M"
celltypes[grep("Cycling G1/S",celltypes)] <- "Neuron associated Cycling G1/S"
celltypes[grep("Neuronal progenitor 4V",celltypes)] <- "Neuronal progenitor 4V"
celltypes[grep("Neuronal progenitor LV",celltypes)] <- "Neuronal progenitor LV"
celltypes[grep("Rspo1",celltypes)] <- "Rspo1+ LV"
celltypes[grep("Progenitor 1",celltypes)] <- "Progenitor 1"
celltypes[grep("Progenitor 2",celltypes)] <- "Progenitor 2"
celltypes[grep("Developing pineal gland (Krt19+)",celltypes)] <- "tmp pineal gland (Krt19+)"
celltypes[grep("Neurons 4V",celltypes)] <- "tmp 4v"
celltypes[grep("Developing pineal gland",celltypes)] <- "Developing pineal gland"
celltypes[grep("Neurons",celltypes)] <- "Neurons"
celltypes[grep("tmp pineal gland (Krt19+)",celltypes)] <- "Developing pineal gland (Krt19+)"
celltypes[grep("tmp 4v",celltypes)] <- "Neurons 4v"
celltypes[grep("Macrophage",celltypes)] <- "Macrophage"
celltypes[grep("Monocyte",celltypes)] <- "Monocyte"
celltypes[grep("Lymphocyte",celltypes)] <- "Lymphocyte"
celltypes[grep("Basophil",celltypes)] <- "Basophil"
celltypes[grep("Neutrophil",celltypes)] <- "Neutrophil"
celltypes[grep("Mast",celltypes)] <- "Mast"
celltypes[grep("DC",celltypes)] <- "DC"
celltypes[grep("CL_0000988\thematopoietic cell\tB",celltypes)] <- "B"

meta_Hchp <- data.frame(cbind(cell=cell_names, celltype=celltypes))
rownames(meta_Hchp) <- cell_names

sce <- SingleCellExperiment(list(counts=counts))
rownames(sce) <- rownames(counts)

meta_Hchp <- meta_Hchp[match(colnames(sce), cell_names),]

colData(sce)  <- DataFrame(meta_Hchp)

sce_filtered  <- sce[calculateAverage(sce)>0.05,]

clusts <- as.numeric(quickCluster(sce_filtered, method = "igraph", min.size = 100))

min.clust <- min(table(clusts))/2
new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce_filtered <- computeSumFactors(sce_filtered, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

sizeFactors(sce) <- sizeFactors(sce_filtered)
#sce <- logNormCounts(sce)
saveRDS(sce, "sce.rds")

###############
library(scran)
library(Seurat)

path2data1 <- "/data1/ivanir/Chp2022/Integrated/"

sce_Hchp <- readRDS(paste0(path2data1, "sce.rds"))

path2data2 <- "/data1/ivanir/Chp2022/ChpMouseCell2021/"

sce_Mchp <- readRDS(paste0(path2data2, "sce.rds"))
rownames(sce_Mchp) <- toupper(rownames(sce_Mchp))

shared_genes <- intersect(rownames(sce_Hchp), rownames(sce_Mchp))

batch <- c(rep("query", ncol(sce_Hchp)), rep("atlas", ncol(sce_Mchp)))

metadata <- data.frame(cbind(cell=c(colnames(sce_Hchp), colnames(sce_Mchp)), batch))
rownames(metadata)   <- c(colnames(sce_Hchp),colnames(sce_Mchp))

### Seurat integration
seurat_integ <- CreateSeuratObject(cbind(counts(sce_Hchp[shared_genes,]), counts(sce_Mchp[shared_genes,])), meta.data = metadata)
seurat_list  <- SplitObject(seurat_integ, split.by = "batch")

## Label transfer
for (i in 1:length(x = seurat_list)) {
   seurat_list[[i]] <- NormalizeData(object = seurat_list[[i]], verbose = TRUE)
   seurat_list[[i]] <- FindVariableFeatures(object = seurat_list[[i]], selection.method = "vst",
        nfeatures = 2000, verbose = TRUE)
}

reference_dataset <- which(names(seurat_list) == "atlas")

seurat_query <- seurat_list[["query"]]
seurat_atlas <- seurat_list[["atlas"]]
seurat_anchors  <- FindTransferAnchors(reference = seurat_atlas, query = seurat_query,
    dims = 1:30)
predictions     <- TransferData(anchorset = seurat_anchors, refdata = colData(sce_Mchp)$celltype,
    dims = 1:30)
seurat_query    <- AddMetaData(object = seurat_query, metadata = predictions)

colData(sce_Hchp)$seurat.prediction <- seurat_query$predicted.id
colData(sce_Hchp)$seurat.max.score  <- seurat_query$prediction.score.max

saveRDS(sce_Hchp, paste0(path2data1, "sce_map2Mchp.rds"))

#############
library(scran)
library(ggplot2)

setwd("/data1/ivanir/Chp2022/Integrated/")
sce <- readRDS("sce_map2Mchp.rds")

celltype_colours <- c(
"Progenitor 1"="#c50068",
"Progenitor 2"="#ff70cc",

"Epithelial cell"="#eda1ff",
"Developing pineal gland"="#6e37bc",

"Neuron associated Cycling G1/S"="#006945",
"Neuron associated Cycling G2/M"="#8f9200",
"Neuronal progenitor 4V"="#54dcb9",
"Neuronal progenitor LV"="#e297a6",
"Neurons"="#4abc34",
"Rspo1+ LV"="#9d241f",

"Lymphocyte"      ="#a96600",

"Mesenchymal cell"="#5fa2ff",

"Below threshold"= "gray"
)

meta <- read.csv("meta_scanpy_updated.csv")
rownames(meta) <- colnames(sce)

df_plot <- data.frame(colData(sce))
umap <- read.csv("umap_layout_scanpy_paga.csv",header=FALSE)
colnames(umap) <- c("UMAP1","UMAP2")
fa <- read.csv("fa2_layout_scanpy_paga.csv",header=FALSE)
colnames(fa) <- c("FA1","FA2")

df_plot <- cbind(df_plot, leiden=meta$leiden, umap, fa)

setwd("/data1/ivanir/Chp2022/ChpMouseCell2021")

reorder_rand <- sample(nrow(df_plot), nrow(df_plot))

ggplot(df_plot[reorder_rand,], aes(x = UMAP1, y = UMAP2, col = factor(seurat.prediction))) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values=celltype_colours[-length(celltype_colours)],name = "Predicted cell type") +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7))) 
ggsave("umap_map2mouse.pdf")

df_plot$seurat_prediction_cutoff <- rep("-",nrow(df_plot))
df_plot$seurat_prediction_cutoff[df_plot$seurat.max.score > .5] <- as.character(df_plot$seurat.prediction[df_plot$seurat.max.score > .5])
plot.index <- order(df_plot$seurat_prediction_cutoff)

ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(seurat_prediction_cutoff))) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values=celltype_colours, name = "Predicted cell type") +
  theme_minimal() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("umap_map2mouse_threshold.pdf")
  
ggplot(df_plot[reorder_rand,], aes(x = FA1, y = FA2, col = factor(seurat.prediction))) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values=celltype_colours,name = "Predicted cell type") +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7))) 
ggsave("fa_map2mouse.pdf")

df_plot$seurat_prediction_cutoff <- rep("-",nrow(df_plot))
df_plot$seurat_prediction_cutoff[df_plot$seurat.max.score > .5] <- as.character(df_plot$seurat.prediction[df_plot$seurat.max.score > .5])
plot.index <- order(df_plot$seurat_prediction_cutoff)

ggplot(df_plot[plot.index,], aes(x = FA1, y = FA2, col = factor(seurat_prediction_cutoff))) +
  geom_point(size = 1, shape = 16) +
  scale_color_manual(values=celltype_colours, name = "Predicted cell type") +
  theme_minimal() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("fa_map2mouse_threshold.pdf")
