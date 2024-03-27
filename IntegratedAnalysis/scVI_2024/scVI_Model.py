import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

import torch
device = torch.device("cpu")

import scvi
import scanpy as sc
import pandas as pd
from skmisc.loess import loess

os.chdir('/data1/ivanir/Chp2022/Integrated/')

adata = sc.read('counts.tab')
adata = adata.T
file = open('cell_names.txt', 'r')
cells  = file.read().splitlines()

file = open('gene_names.txt', 'r')
genes   = file.read().splitlines()

adata.obs = pd.read_csv('meta.tab', sep ='\t', low_memory=False)
adata.obs_names = cells
adata.var_names = genes

adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # freeze the state in `.raw`

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="batch",
)

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["batch", "platform"],
    #continuous_covariate_keys=["day"],
)
model = scvi.model.SCVI(adata)
model

model.train()

model_dir = os.path.join('/data1/ivanir/Chp2022/Integrated/', "scvi_model")
model.save(model_dir, overwrite=True)

#model = scvi.model.SCVI.load(model_dir, adata=adata)

SCVI_LATENT_KEY = "X_scVI"

latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
latent.shape

sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=0.3)

sc.pl.umap(
    adata,
    color=["platform"],
    frameon=False,
)

sc.write("scvi_adata.h5ad",adata)

adata = sc.read("scvi_adata.h5ad")

sc.tl.draw_graph(adata)

sc.pl.draw_graph(
    adata,
    color=["platform"],
    frameon=False,
)

sc.write("scvi_adata.h5ad",adata)

