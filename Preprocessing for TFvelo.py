import scvelo as sv
import scanpy as sc
import anndata as ad
import muon as mu
import mudata as md
import numpy as np
from sklearn.impute import SimpleImputer



adata = sc.read_h5ad("./Desktop/Human_DCA/anndata/DCA_Mac.h5ad")

adata.obs["seurat_clusters"] = adata.obs["seurat_clusters"].astype(str)

sc.pl.umap(adata, color="seurat_clusters")

adata.obs["clusters"] = adata.obs["seurat_clusters"]
adata.uns["clusters_colors"] = adata.uns["seurat_clusters_colors"]

#Mac
new_colors = np.array(adata.uns['clusters_colors'])
new_colors[[0, 1, 2]] = np.array(['#f8766d', '#b79f00', '#00ba38'])  # outliers / grey
adata.uns['clusters_colors'] = new_colors
sc.pl.umap(adata, color="clusters")


imputer = SimpleImputer(strategy='constant',fill_value = 0.0)
adata.X = imputer.fit_transform(adata.X)


#sc.tl.leiden(adata)
sc.pl.umap(adata, color="seurat_clusters")

adata_loom1 = sv.read('./Desktop/Human_DCA/CellRanger210615/1SA/velocyto1/1SA.loom', cache=True)
adata_loom2 = sv.read('./Desktop/Human_DCA/CellRanger210615/2ACS/velocyto/2ACS.loom', cache=True) 


sv.utils.clean_obs_names(adata)
sv.utils.clean_obs_names(adata_loom1)
sv.utils.clean_obs_names(adata_loom2)
adata1 = sv.utils.merge(adata, adata_loom1)
adata2 = sv.utils.merge(adata, adata_loom2)

#adata_test =sv.utils.merge(adata,adata_loom2,adata_loom1)

#All
common_genes = set(adata1.var_names)
common_genes.intersection_update(adata2.var_names)

adata1 = adata1[:, list(common_genes)]
adata2 = adata2[:, list(common_genes)]

var1 = adata1.var.copy()
var2 = adata2.var.copy()
#var3 = adata3.var.copy()
#var4 = adata4.var.copy()

for col in var1.columns:
    if col not in var2.columns:
        var2[col] = np.nan

for col in var2.columns:
    if col not in var1.columns:
        var1[col] = np.nan

for col in var3.columns:
    if col not in var1.columns:
        var1[col] = np.nan
    if col not in var2.columns:
        var2[col] = np.nan
    if col not in var4.columns:
        var4[col] = np.nan

for col in var4.columns:
    if col not in var1.columns:
        var1[col] = np.nan
    if col not in var2.columns:
        var2[col] = np.nan
    if col not in var3.columns:
        var3[col] = np.nan

var_combined = var1.combine_first(var2)
#var_combined = var1.combine_first(var2).combine_first(var3).combine_first(var4)

adatas = {"1SA": adata1, "2ACS":adata2 }
adatas = ad.concat(adatas, label="dataset_name",  join="outer")

adatas.var =  var_combined

sc.pl.umap(adatas, color="clusters")

adatas.obs["clusters"] = adatas.obs["seurat_clusters"]
adatas.uns["clusters_colors"] = adatas.uns["seurat_clusters_colors"]

#Mye
new_colors = np.array(adatas.uns['clusters_colors'])
new_colors[[0, 1, 2, 3, 4, 5]] = np.array(['#f8766d', '#b79f00', '#00ba38', '#00bfc4', '#619cff', '#f564e3'])  # outliers / grey
adatas.uns['clusters_colors'] = new_colors
sc.pl.umap(adatas, color="clusters")

#Mac
new_colors = np.array(adatas.uns['clusters_colors'])
new_colors[[0, 1, 2]] = np.array(['#f8766d', '#b79f00', '#00ba38'])  # outliers / grey
adatas.uns['clusters_colors'] = new_colors
sc.pl.umap(adatas, color="clusters")

adatas.X = imputer.fit_transform(adatas.X)

adatas.write("./Desktop/Human_DCA/anndata/DCA_all_Myeloid_with_loom_240711.h5ad", compression='gzip')
adatas.write("./Desktop/TFvelo/data/DCA/DCA_all_Myeloid_with_loom_240711.h5ad", compression='gzip')


adatas_sa = adatas[adatas.obs.stim == "sa"]
imputer = SimpleImputer(strategy='constant',fill_value = 0.0)
adatas_sa.X = imputer.fit_transform(adatas_sa.X)

sc.pl.umap(adatas_sa, color="clusters")

adatas_uapmi = adatas[adatas.obs.stim == "uapmi"]
imputer = SimpleImputer(strategy='constant',fill_value = 0.0)
adatas_uapmi.X = imputer.fit_transform(adatas_uapmi.X)

adatas_sa.write("./Desktop/Human_DCA/anndata/DCA_sa_Myeloid_with_loom_240918.h5ad", compression='gzip')
adatas_sa.write("./Desktop/TFvelo/data/DCA/DCA_sa_Myeloid_with_loom_240918.h5ad", compression='gzip')

adatas_uapmi.write("./Desktop/Human_DCA/anndata/DCA_uapmi_Myeloid_with_loom_240918.h5ad", compression='gzip')
adatas_uapmi.write("./Desktop/TFvelo/data/DCA/DCA_uapmi_Myeloid_with_loom_240918.h5ad", compression='gzip')

#without loom
adatas_sa = adata[adata.obs.stim == "sa"]
imputer = SimpleImputer(strategy='constant',fill_value = 0.0)
adatas_sa.X = imputer.fit_transform(adatas_sa.X)

sc.pl.umap(adatas_sa, color="clusters")

adatas_uapmi = adata[adata.obs.stim == "uapmi"]
imputer = SimpleImputer(strategy='constant',fill_value = 0.0)
adatas_uapmi.X = imputer.fit_transform(adatas_uapmi.X)

adatas_sa.var.rename(columns={"_index": "index"}, inplace=True)
adatas_uapmi.var.rename(columns={"_index": "index"}, inplace=True)

if adatas_uapmi.raw is not None:
    adatas_uapmi.raw.var.rename(columns={"_index": "index"}, inplace=True)
if adatas_sa.raw is not None:
    adatas_sa.raw.var.rename(columns={"_index": "index"}, inplace=True)

adatas_sa.write("./Desktop/Human_DCA/anndata/DCA_sa_Mac_240918.h5ad", compression='gzip')
adatas_sa.write("./Desktop/TFvelo/data/DCA/DCA_sa_Mac_240918.h5ad", compression='gzip')

adatas_uapmi.write("./Desktop/Human_DCA/anndata/DCA_uapmi_Mac_240918.h5ad", compression='gzip')
adatas_uapmi.write("./Desktop/TFvelo/data/DCA/DCA_uapmi_Mac_240918.h5ad", compression='gzip')

