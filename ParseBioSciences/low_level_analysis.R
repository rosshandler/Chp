library(scran)
library(irlba)
library(Rtsne)
library(Matrix)
library(ggplot2)
library(biomaRt)
library(viridis)
library(scDblFinder)

library(umap)
library(reticulate)
use_condaenv(condaenv="scanpy-p3.9")

umap = import('umap')

path2data  <- '/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/Combined/all-well/DGE_filtered/'

counts    <- t(readMM(paste0(path2data, "DGE.mtx")))
genes     <- read.csv(paste0(path2data, "all_genes.csv"))
metadata  <- read.csv(paste0(path2data, "cell_metadata.csv"))

sample_info <- read.table(paste0(path2data,'sample_info.tab'), header = TRUE)

lib.sizes <- colSums(counts)
ngenes    <- colSums(counts > 0)

seq_library <- rep(NA, nrow(metadata))
sample_bc1_well <- rep(NA, nrow(metadata))        
sample_number   <- rep(NA, nrow(metadata))
sample_name     <- rep(NA, nrow(metadata))

samples <- unique(sample_info$sample_bc1_well)
for (i in 1:length(samples)){
  sample_bc1_well[metadata$bc1_well %in% samples[i]] <- sample_info$sample_bc1_well[i]
  sample_number[metadata$bc1_well %in% samples[i]]   <- sample_info$sample_number[i]
  sample_name[metadata$bc1_well %in% samples[i]]     <- sample_info$sample_name[i]
}
seq_library[grep("s1",metadata$bc_wells)] <- "DNAA007"
seq_library[grep("s2",metadata$bc_wells)] <- "DNAA008"

day <- gsub("n.*","",sample_name)
metadata <- data.frame(cbind(metadata,lib.sizes, cbind(sample_bc1_well,sample_number,sample_name,seq_library, day)))

plot_df <- metadata

pdf("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/UMIsByReplicate.pdf")
ggplot(plot_df, aes (x = factor(sample_replicate), y = as.numeric(lib.sizes))) +
  geom_boxplot() +
  theme_bw() + 
  labs(x = "Batch", y = "Number of UMIs") +
  scale_y_log10(breaks = c(1000, 5000, 10000, 50000, 100000),
    labels = c("1,000", "5,000", "10,000", "50,000", "100,000"))
dev.off()

pdf("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/cell_complexity.pdf")
qplot(lib.sizes, ngenes, col = ifelse(ngenes < 500, "drop", "keep")) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  labs(x = "UMI count", y = "Number of expressed genes") +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
dev.off()

ensembl <- useEnsembl(biomart = "ensembl",  dataset = "hsapiens_gene_ensembl",mirror="useast")

gene_map  <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
  filters = "hgnc_symbol", values = genes$gene_name, mart = ensembl)

mt.index    <- gene_map$chromosome_name == "MT"
mt.counts   <- counts[which(genes$gene_name %in% gene_map$hgnc_symbol[mt.index]), ]
mt.fraction <- colSums(mt.counts)/lib.sizes

mt.p   <- pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
mt.lim <- min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 0.05)])

#Threhdold
mt.lim
#[1] 0.08991269

mt.lim <- min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 0.001)])

#Threhdold
mt.lim
#[1] 0.1183915

pdf("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/mtreadfraction.pdf")
qplot(lib.sizes, mt.fraction, col = ifelse(mt.fraction>mt.lim, "drop", "keep")) +
  scale_x_log10() +
  labs(x = "UMI count", y = "MT read fraction") +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
dev.off()

filter.cell <- mt.fraction > mt.lim
plot_df <- data.frame(lib.size = lib.sizes, mt.fraction, filter.cell, sample = samples)

ggplot(plot_df, aes(fill=filter.cell, y=1, x=sample)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E", name = "Filter cell") + 
    ylab("Fraction of cells removed by stress signals") + xlab("Sample")
ggsave("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/mt_fraction_reomval.pdf")

mt.lim =.2

pdf("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/mtreadfraction_manual.pdf")
qplot(lib.sizes, mt.fraction, col = ifelse(mt.fraction>mt.lim, "drop", "keep")) +
  scale_x_log10() +
  labs(x = "UMI count", y = "MT read fraction") +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
dev.off()

filter.cell <- mt.fraction > mt.lim
plot_df <- data.frame(lib.size = lib.sizes, mt.fraction, filter.cell, sample = samples)

ggplot(plot_df, aes(fill=filter.cell, y=1, x=sample)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E", name = "Filter cell") + 
    ylab("Fraction of cells removed by stress signals") + xlab("Sample")
ggsave("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/mt_fraction_reomval_manual.pdf")

counts      <- counts[, mt.fraction < mt.lim]
metadata    <- metadata[mt.fraction < mt.lim,]
sample_bc1_well <- sample_bc1_well[mt.fraction < mt.lim]

mt.fraction <- mt.fraction[mt.fraction < mt.lim]
lib.sizes <- colSums(counts)
n.genes   <- colSums(counts>0)

plot_df <- data.frame(lib.size = lib.sizes, sample = sample_bc1_well)
level_order <- samples
plot_df$sample <- factor(plot_df$sample, levels=level_order)

ggplot(plot_df, aes(x = factor(sample), y= lib.size)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Sample", y = "Number of UMIs")
ggsave("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/UMIsBySample.pdf")

cell     <- paste0("cell_",1:nrow(metadata))
metadata <- cbind(cell, metadata, mt.fraction)

sce <- SingleCellExperiment(list(counts=counts),colData=DataFrame(metadata))
rownames(sce) <- genes$gene_id

rownames(genes) <- rownames(sce)
rowData(sce) <- DataFrame(genes)

colnames(sce) <- cell
rownames(metadata) <- colnames(sce)
colData(sce) <- DataFrame(metadata)

write.table(metadata, file = paste0(path2data, "cell_metadata.tab"), row.names = F, col.names = T, quote = F, sep = "\t")

lib.sizes <- colSums(counts(sce))
sce_filt  <- sce[calculateAverage(sce)>0.05,]

clusts <- as.numeric(quickCluster(sce_filt, method = "igraph", min.size = 100))

min.clust <- min(table(clusts))/2
new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce_filt <- computeSumFactors(sce_filt, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

sizeFactors(sce) <- sizeFactors(sce_filt)

ggplot(data = data.frame(X = lib.sizes, Y = sizeFactors(sce)), mapping = aes(x = X, y = Y)) +
  geom_point() +
  scale_x_log10(breaks = c(500, 2000, 5000, 10000, 30000), labels = c("5,00", "2,000", "5,000", "10,000", "30,000") ) +
  scale_y_log10(breaks = c(0.2, 1, 5)) +
  theme_minimal() +
  theme(text = element_text(size=20))  +
  labs(x = "Number of UMIs", y = "Size Factor")
ggsave("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/sizefactors.pdf")

write.table(sizeFactors(sce_filt), quote = F, col.names = F, row.names = F,
  file = paste0(path2data, "sizefactors.tab"))

library(BiocParallel)

bp <- MulticoreParam(12, RNGseed=1234)
bpstart(bp)
sce <- scDblFinder(sce, samples="sample_number", dbr=.03, dims=30, BPPARAM=bp)
bpstop(bp)
table(sce$scDblFinder.class)
#singlet doublet 
#  13059     540  

saveRDS(sce,paste0(path2data,"sce.rds"))

lib.size <- colSums(counts(sce))

df_plot <- data.frame(
 lib.size = log10(lib.size),
 doublet  = colData(sce)$scDblFinder.class
)

ggplot(df_plot, aes(x=doublet, y=lib.size)) + 
  geom_boxplot() +
  labs(x = "", y = "log10(Library size)") 
ggsave("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/library_size_doublets.pdf")

sce <- sce[calculateAverage(sce)>0.05,]

sce <- logNormCounts(sce)

decomp  <- modelGeneVar(sce)
hvgs    <- rownames(decomp)[decomp$FDR < 0.1]
pca     <- prcomp_irlba(t(logcounts(sce[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce)
tsne <- Rtsne(pca$x, pca = FALSE, check_duplicates = FALSE)
layout  <- umap(pca$x, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs", "min_dist"), n_neighbors=15, min_dist=.99)

df_plot <- data.frame(
 metadata,
 tSNE1    = tsne$Y[, 1],
 tSNE2    = tsne$Y[, 2], 
 doublet = colData(sce)$scDblFinder.class,
 UMAP1 = layout$layout[,1],
 UMAP2 = layout$layout[,2] 
)

plot.index <- order(df_plot$doublet)
ggplot(df_plot[plot.index,], aes(x = tSNE1, y = tSNE2, col = factor(doublet))) +
  geom_point(size = 0.4) +
  scale_color_manual(values=c("gray","#0169c1"), name = "") +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/tsne_doublets.pdf")

samples <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12")
days <- c("14",  "20", "25", "70", "95", "130")

level_order <- days
df_plot$day <- factor(df_plot$day,levels=level_order)

level_order <- samples
df_plot$sample <- factor(df_plot$bc1_well, levels=level_order)

plot.index <- sample(nrow(df_plot))

ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = sample)) +
  geom_point(size = 1) +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/umap_by_library.pdf")

ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = day)) +
  geom_point(size = 1) +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/umap_by_day.pdf")

sample_colours <- c(
"A1"="#c85d3e",
"A2"="#46afc8",
"A3"="#d84165",
"A4"="#55b871",
"A5"="#ca4ba8",
"A6"="#9ab341",
"A7"="#8e5fcc",
"A8"="#617b39",
"A9"="#6d7cc7",
"A10"="#c69042",
"A11"="#d78ac0",
"A12"="#a85171")

ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = bc1_well)) +
  geom_point(size = 1) +
  scale_color_manual(values=sample_colours, name = "Sample") +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("/data1/ivanir/Chp2022/ParseBS/newvolume/analysis/QC/umap_by_sample.pdf")
  
