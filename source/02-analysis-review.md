---
jupyter:
  celltoolbar: Slideshow
  jupytext:
    cell_metadata_json: true
    formats: ipynb,md,py:percent
    notebook_metadata_filter: all
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
  language_info:
    codemirror_mode:
      name: ipython
      version: 3
    file_extension: .py
    mimetype: text/x-python
    name: python
    nbconvert_exporter: python
    pygments_lexer: ipython3
    version: 3.10.7
  rise:
    scroll: true
    theme: black
  toc-autonumbering: true
  toc-showcode: false
  toc-showmarkdowntxt: false
  widgets:
    application/vnd.jupyter.widget-state+json:
      state: {}
      version_major: 2
      version_minor: 0
---

<!-- #region {"slideshow": {"slide_type": "slide"}, "tags": []} -->
# Analysis review
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
## Overview
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
This notebook follows the tutorial by [mousepixels/sanbomics](https://github.com/mousepixels/sanbomics/blob/main/single_cell_analysis_complete_class.ipynb), which has an accompanying [screencast](https://youtu.be/uvyG9yLuNSE?t=319).

Analysis is illustrated with single-nucleus RNA sequencing data from the following paper <cite data-cite="Melms2021-bj">Melms et al. (2021)</cite>

> Melms JC, Biermann J, Huang H, Wang Y, Nair A, Tagore S, et al.
A molecular single-cell lung atlas of lethal COVID-19.
Nature. 2021;595: 114–119. [doi:10.1038/s41586-021-03569-0](https://doi.org/10.1038/s41586-021-03569-0)

This paper examined 116,000 nuclei from the lungs of nineteen patients who underwent autopsy following death in association with COVID-19. Findings reported in the abstract of the paper include:

1. activated monocyte-derived macrophages and alveolar macrophages
1. impaired T cell activation
1. monocyte/macrophage-derived interleukin-1β and epithelial cell-derived interleukin-6
1. alveolar type 2 cells adopted an inflammation-associated transient progenitor cell state and failed to undergo full transition into alveolar type 1 cells
1. expansion of CTHRC1+ pathological fibroblasts
1. protein activity and ligand–receptor interactions suggest putative drug targets

This notebook makes extensive use of <cite data-cite="Wolf2018-nu">Wolf et al. (2018)</cite> and <cite data-cite="Lopez2018-em">Lopez et al. (2018)</cite> including updates that have been made to the underlying software packages, [scanpy](https://github.com/scverse/scanpy) and [scvi-tools](https://github.com/scverse/scvi-tools), since their initial publication.
<!-- #endregion -->

<!-- #region {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"slide\"} jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "slide"}, "tags": []} -->
## Setup
<!-- #endregion -->

<!-- #region {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"subslide\"} jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Import libraries
<!-- #endregion -->

```python tags=[] slideshow={"slide_type": "fragment"} jupyter={"outputs_hidden": true}
from inspect import getmembers
from pprint import pprint
from types import FunctionType

import numpy as np
import pickle
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
```

<!-- #region {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"subslide\"} jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Setup plotting
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
import matplotlib.font_manager
import matplotlib.pyplot as plt

# import matplotlib_inline
```

```python slideshow={"slide_type": "fragment"} tags=[]
# fonts_path = "/usr/share/texmf/fonts/opentype/public/lm/" #ubuntu
# fonts_path = "~/Library/Fonts/" # macos
fonts_path = "/usr/share/fonts/OTF/"  # arch
# user_path = "$HOME/" # user
# fonts_path = user_path + "fonts/latinmodern/opentype/public/lm/"  # home
matplotlib.font_manager.fontManager.addfont(fonts_path + "lmsans10-regular.otf")
matplotlib.font_manager.fontManager.addfont(fonts_path + "lmroman10-regular.otf")
```

```python slideshow={"slide_type": "fragment"} tags=[]
# https://stackoverflow.com/a/36622238/446907
%config InlineBackend.figure_formats = ['svg']
```

```python slideshow={"slide_type": "fragment"} tags=[]
plt.style.use("default")  # reset default parameters
# https://stackoverflow.com/a/3900167/446907
plt.rcParams.update(
    {
        "font.size": 16,
        "font.family": ["sans-serif"],
        "font.serif": ["Latin Modern Roman"] + plt.rcParams["font.serif"],
        "font.sans-serif": ["Latin Modern Sans"] + plt.rcParams["font.sans-serif"],
    }
)
```

<!-- #region {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"subslide\"} jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Utility functions
<!-- #endregion -->

```python tags=[] slideshow={"slide_type": "fragment"}
def attributes(obj):
    """
    get object attributes
    """
    disallowed_names = {
        name for name, value in getmembers(type(obj)) if isinstance(value, FunctionType)
    }
    return {
        name: getattr(obj, name)
        for name in dir(obj)
        if name[0] != "_" and name not in disallowed_names and hasattr(obj, name)
    }


def print_attributes(obj):
    """
    print object attributes
    """
    pprint(attributes(obj))
```

<!-- #region {"slideshow": {"slide_type": "slide"}, "tags": []} -->
## Import data
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Here we review how the data were downloaded, and proceed to import and inspect the data.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Data download
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Data with GEO accession [GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524) was downloaded using [./data/download_geo_data.sh](./data/download_geo_data.sh) with parameters

```bash
./download_geo_data.sh \
       -a GSE132771 \
       -f 'ftp.*RAW.*' \
       -j '..|.supplementary_files?|..|.url?|select(length>0)'
```

A skeleton of this script that may work in this case is

```bash
!/usr/bin/env bash

#-- debugging (comment to reduce stderr output)
#-- https://wiki.bash-hackers.org/scripting/debuggingtips
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
set -o xtrace

# get metadata
# Melms JC, Biermann J, Huang H, Wang Y, Nair A, Tagore S, et al.
# A molecular single-cell lung atlas of lethal COVID-19.
# Nature. 2021;595: 114–119. doi:10.1038/s41586-021-03569-0
# GSE171524
ffq -l 1 -o GSE171524.json GSE171524

# download raw data
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171524/suppl/GSE171524_RAW.tar

# list contents
tar -tvf GSE171524_RAW.tar

# untar
mkdir -p GSE171524 && \
tar -xvf GSE171524_RAW.tar -C GSE171524
```
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Data load
<!-- #endregion -->

See the [documentation for scanpy read csv](https://scanpy.readthedocs.io/en/latest/generated/scanpy.read_csv.html) which returns an [AnnData object](https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.html#anndata.AnnData).

```python slideshow={"slide_type": "fragment"} tags=[]
adata = None
adata = sc.read_csv(
    "data/GSE171524/supplementary/GSM5226574_C51ctr_raw_counts.csv.gz"
).T
adata
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Note the `scanpy.read_csv` function accepts gzipped files.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Data properties
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
type(adata)
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.X.shape
```

```python tags=[] slideshow={"slide_type": "subslide"}
print_attributes(adata)
```

```python tags=[] slideshow={"slide_type": "subslide"}
adata.obs
```

Gene names are saved as `adata.var`.

```python tags=[] slideshow={"slide_type": "fragment"}
adata.var
```

```python tags=[] slideshow={"slide_type": "subslide"}
adata.obs_names
```

```python tags=[] slideshow={"slide_type": "fragment"}
adata.var_names
```

There are no layers in this data set.

```python tags=[] slideshow={"slide_type": "subslide"}
adata.layers
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
There are no multidimensional observations or variables.
<!-- #endregion -->

```python tags=[] slideshow={"slide_type": "fragment"}
print(adata.obsm)
print(adata.varm)
print(adata.obsp)
print(adata.varp)
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
The data appears to contain reads mapped to 6099 cell-associated barcodes and 34546 RNA molecule-associated features.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "slide"}, "tags": []} -->
## Doublet removal
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Filter transcripts by minimum number of cells with non-zero counts
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
We may choose to filter out transcript types that are detected in a relatively small number of cells. The optimum threshold is not known. Here we use 10 as a base case.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
See the [documentation for scanpy pre-processing filter-genes](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.filter_genes.html).
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
sc.pp.filter_genes(adata, min_cells=10)
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
Following filtration there are 19896 transcript types remaining that are present in at least 10 cells. This means 14650 transcript types were removed at the threshold of 10. Notice that an annotation named `n_cells` has been added to the genes to indicate the number of cells with non-zero values for that transcript type.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata.var
```

```python
adata.var.describe()
```

```python tags=[] slideshow={"slide_type": "subslide"}
sns.histplot(adata.var['n_cells'])
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Select highly variable genes
<!-- #endregion -->

It is common to select transcript types with high variability among the cell population under the assumption that this will help to focus on features that distinguish the cells. Again there is no perfect threshold. Here we select the 2000 highest variability genes.

```python tags=[]
adata
```

See the [documentation for scanpy pre-processing highly-variable-genes](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.highly_variable_genes.html).

```python tags=[]
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor="seurat_v3")
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
We have filtered down to 2000 transcript types and added 5 annotations to the genes including a binary variable indicating membership in the highly-variable class, the ranking among highly variable genes, and the mean, variance, and normalized variance for each gene across cells.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata.var
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Doublet removal
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Here we `scvi-tools` model as [an interface](https://docs.scvi-tools.org/en/stable/api/reference/scvi.external.SOLO.html) to [solo](https://github.com/calico/solo) from <cite data-cite="Bernstein2020-am">Bernstein et al. (2020)</cite>. `scvi` adds some unstructured data to interface with `torch`.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata.uns
```

```python tags=[] slideshow={"slide_type": "fragment"}
scvi.model.SCVI.setup_anndata(adata)
```

```python tags=[] slideshow={"slide_type": "fragment"}
adata.uns
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
Train the `scvi` model.
<!-- #endregion -->

```python tags=[] slideshow={"slide_type": "fragment"}
vae = scvi.model.SCVI(adata)
vae.train()
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
After training the model `adata.obs` now contains `_scvi_batch` and `scvi_labels` annotations.
<!-- #endregion -->

```python tags=[] slideshow={"slide_type": "fragment"}
adata
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.obs
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
Train the `solo` model.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
solo = scvi.external.SOLO.from_scvi_model(vae)
solo.train()
```

```python slideshow={"slide_type": "fragment"} tags=[]
type(solo)
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
We can now use the `solo` model to add a prediction annotation for doublet or singlet to our data.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
df = solo.predict()
df["prediction"] = solo.predict(soft=False)

df.index = df.index.map(lambda x: x[:-2])

df
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
Counting the predicted doublets and singlets reveals about 20% doublets.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
df.groupby("prediction").count()
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
We can assess the magnitude of the prediction by taking the difference between the doublet and singlet scores.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
df["dif"] = df.doublet - df.singlet
df
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
Plotting the distribution of doublet-singlet score differences we see there are a large fraction that marginally favor doublet to singlet.
<!-- #endregion -->

```python tags=[] slideshow={"slide_type": "fragment"}
sns.histplot(df[df.prediction == "doublet"], x="dif")
```

Since we would like to avoid unnecessarily discarding barcodes, we will arbitrarily retain those with a score from zero to one (however this is not intended to be principled).

```python
doublets = df[(df.prediction == 'doublet') & (df.dif > 1)]
doublets
```

```python tags=[]
sns.histplot(doublets[doublets.prediction == "doublet"], x="dif")
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Save/load checkpoint 
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Save current variables to file.
<!-- #endregion -->

```python tags=[] slideshow={"slide_type": "fragment"}
pickle.dump([adata, df, doublets, solo, vae], open("session_09282022.p", "wb"))
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Variables can be reloaded if necessary.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata, df, doublets, solo, vae = pickle.load(open("session_09282022.p", "rb"))
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Reload data and filter doublets
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata = sc.read_csv(
    "data/GSE171524/supplementary/GSM5226574_C51ctr_raw_counts.csv.gz"
).T
adata
```

Annotate `adata` with doublet column.

```python slideshow={"slide_type": "fragment"} tags=[]
adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.obs
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
Use doublet column to filter doublets.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata = adata[~adata.obs.doublet]
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
We filtered about 7% of the barcodes as putative doublets.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "slide"}, "tags": []} -->
## Preprocessing
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Filter mitochondrial transcripts
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
We can add a binary annotation to indicate presence of mitochondrial transcripts.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata.var['mt'] = adata.var.index.str.startswith('MT-')
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.var
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Filter ribosomal transcripts
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
We can download a list of ribosomal genes from msigdb.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
ribo_genes = pd.read_table(ribo_url, skiprows=2, header = None)
ribo_genes
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
We add a binary annotation to indicate presence in the list of ribosomal genes.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.var
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Calculate QC metrics
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Recall the observation/barcode annotation currently just indicates doublets based on our solo model.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata.obs
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
We can calculate standard quality control metrics with `scanpy`. See the [documentation for scanpy calculate_qc_metrics](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.calculate_qc_metrics.html). These include the number of 
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
Sorting genes by the number of cells with non-zero counts, we see a number of genes that were not found in association with any barcode.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata.var.sort_values("n_cells_by_counts")
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
Sorting barcodes by total counts we see a minimum around 400 suggesting a previously applied filter upstream of this pre-processing.
<!-- #endregion -->

```python tags=[] slideshow={"slide_type": "fragment"}
adata.obs.sort_values("total_counts")
```

```python slideshow={"slide_type": "subslide"} tags=[]
sns.jointplot(
    data=adata.obs,
    x="total_counts",
    y="n_genes_by_counts",
    kind="hex",
)
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Filter genes by counts
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
We can filter genes found in some minimum number of cells, in this case we arbitrarily choose three.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
sc.pp.filter_genes(adata, min_cells=3)
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.var.sort_values("n_cells_by_counts")
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Filter cells by number of genes
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
As a matter of completeness we filter cells with fewer than 200 genes; however, we see the minimum number of genes in a cell is already greater than 200.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata.obs.sort_values('n_genes_by_counts')
```

```python slideshow={"slide_type": "subslide"} tags=[]
sc.pp.filter_cells(adata, min_genes=200)
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Review and filter QC outliers
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
We can plot distributions of the number of genes by counts and total counts as well as the percentage of mitochondrial and ribosomal genes.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True)
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
If we choose to filter cells beyond the 98th percentile, we can compute the bound via the numpy quantile function.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
upper_lim
```

```python tags=[] slideshow={"slide_type": "fragment"}
adata = adata[adata.obs.n_genes_by_counts < upper_lim]
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
This filters a few more than 100 additional cells.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata.obs
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
Finally, we filter mitochondrial (at 20%) and ribosomal (at 2%) outliers.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata = adata[adata.obs.pct_counts_mt < 20]
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata = adata[adata.obs.pct_counts_ribo < 2]
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
We end our preprocessing with 5518 cells and 24136 genes.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Save/load checkpoint 
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Save current variables to file.
<!-- #endregion -->

```python tags=[] slideshow={"slide_type": "fragment"}
pickle.dump([adata, df, doublets, solo, vae], open("adata_post_processed.p", "wb"))
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Variables can be reloaded if necessary.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata, df, doublets, solo, vae = pickle.load(open("adata_post_processed.p", "rb"))
```

## Normalization

```python
adata.X.sum(axis=1)
```

```python
sc.pp.normalize_total(adata, target_sum=1e4) 
```

```python
adata.X.sum(axis=1)
```

```python
sc.pp.log1p(adata)
```

```python
adata.X.sum(axis=1)
```

```python tags=[]
adata.raw = adata
```
