###### conda activate scvelo
import os
import numpy as np
import loompy as lp
import pandas as pd
import scanpy as sc

import scipy.sparse as sparse
import scipy.io as sio
import scipy.stats as stats

os.chdir('//data1/ivanir/Chp2022/10xGen/analysis/')

lp.combine(
["/data2/mlancast/cruk/newfastq/ChoroidPlexusC/velocyto/ChoroidPlexusC.loom",
 "/data2/mlancast/cruk/newfastq/ChoroidPlexusD/velocyto/ChoroidPlexusD.loom",
 "/data2/mlancast/cruk/newfastq/ChoroidPlexusE/velocyto/ChoroidPlexusE.loom",
 "/data2/mlancast/cruk/SLX-20683/fastqs/ChoroidPlexusH1D19/velocyto/ChoroidPlexusH1D19.loom",
 "/data2/mlancast/cruk/SLX-20683/fastqs/ChoroidPlexusH1D33/velocyto/ChoroidPlexusH1D33.loom",
 "/data2/mlancast/cruk/SLX-20683/fastqs/ChoroidPlexusH1D80/velocyto/ChoroidPlexusH1D80.loom",
 "/data2/mlancast/cruk/SLX-20683/fastqs/ChoroidPlexusH9D63/velocyto/ChoroidPlexusH9D63.loom",
 "/data1/ivanir/Chp2022/10xGen/fastqs/ChPH1D19_A/velocyto/ChPH1D19_A.loom",
 "/data1/ivanir/Chp2022/10xGen/fastqs/ChPH1D19_B/velocyto/ChPH1D19_B.loom"],
 "ChoroidPlexus10X.loom")

adata = sc.read('ChoroidPlexus10X.loom', sparse=True)
adata.var_names_make_unique()

sio.mmwrite("unspliced_matrix.mtx",adata.layers['unspliced'].T, field='integer')
sio.mmwrite("spliced_matrix.mtx",adata.layers['spliced'].T, field='integer')

np.savetxt('barcodes_loom.csv', adata.obs_names.values,fmt='%5s', delimiter=',')
np.savetxt('gene_ids_loom.csv', adata.var_names.values,fmt='%5s', delimiter=',')
