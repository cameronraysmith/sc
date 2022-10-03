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

```python slideshow={"slide_type": "fragment"} tags=[]
from inspect import getmembers
from pprint import pprint
from types import FunctionType
from datetime import datetime
import os

import numpy as np
import pickle
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
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

```python slideshow={"slide_type": "fragment"} tags=[]
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

### Set output folder for session

```python tags=[]
now = datetime.now()
timestamp = now.strftime("%Y%m%d_%H%M%S")
timestamp = "20220930_215746"
```

```python tags=[]
output_directory = f"output/{timestamp}"
```

```python tags=[]
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
    print(f"created {output_directory}")
```

```python tags=[]
print(f"using {output_directory}")
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

```python slideshow={"slide_type": "subslide"} tags=[]
print_attributes(adata)
```

```python slideshow={"slide_type": "subslide"} tags=[]
adata.obs
```

Gene names are saved as `adata.var`.

```python slideshow={"slide_type": "fragment"} tags=[]
adata.var
```

```python slideshow={"slide_type": "subslide"} tags=[]
adata.obs_names
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.var_names
```

There are no layers in this data set.

```python slideshow={"slide_type": "subslide"} tags=[]
adata.layers
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
There are no multidimensional observations or variables.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
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

```python slideshow={"slide_type": "subslide"} tags=[]
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

```python slideshow={"slide_type": "fragment"} tags=[]
scvi.model.SCVI.setup_anndata(adata)
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.uns
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
Train the `scvi` model.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
vae = scvi.model.SCVI(adata)
vae.train()
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
After training the model `adata.obs` now contains `_scvi_batch` and `scvi_labels` annotations.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
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

```python slideshow={"slide_type": "fragment"} tags=[]
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

```python slideshow={"slide_type": "fragment"} tags=[]
pickle.dump([df, doublets, solo, vae], open(f"{output_directory}/auxiliary.p", "wb"))
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Variables can be reloaded if necessary.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
df, doublets, solo, vae = pickle.load(open(f"{output_directory}/auxiliary.p", "rb"))
```

```python tags=[]
adata.write(f"{output_directory}/adata_doublets.h5ad", compression="gzip")
```

```python tags=[]
bdata_doublets = None
bdata_doublets = sc.read_h5ad(f"{output_directory}/adata_doublets.h5ad")
bdata_doublets
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

```python slideshow={"slide_type": "fragment"} tags=[]
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

```python slideshow={"slide_type": "fragment"} tags=[]
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
Save current version of post processed data to file.
<!-- #endregion -->

```python tags=[]
adata.write(f"{output_directory}/adata_post_processed.h5ad", compression="gzip")
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Variables can be reloaded if necessary.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
df, doublets, solo, vae = pickle.load(open(f"{output_directory}/auxiliary.p", "rb"))
```

```python
bdata_post_processed = None
bdata_post_processed = sc.read_h5ad(f"{output_directory}/adata_post_processed.h5ad")
bdata_post_processed
```

<!-- #region {"slideshow": {"slide_type": "slide"}, "tags": []} -->
## Normalization
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata.X.sum(axis=1)
```

```python slideshow={"slide_type": "fragment"} tags=[]
sc.pp.normalize_total(adata, target_sum=1e4) 
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.X.sum(axis=1)
```

```python slideshow={"slide_type": "fragment"} tags=[]
sc.pp.log1p(adata)
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.X.sum(axis=1)
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.raw = adata
```

<!-- #region {"slideshow": {"slide_type": "slide"}, "tags": []} -->
## Clustering
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Estimate expression variability
<!-- #endregion -->

```python slideshow={"slide_type": "skip"} tags=[]
adata.var
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Here we use the `scanpy` pre-processing [function for highly-variable genes](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.highly_variable_genes.html) to annotate genes by their expression variability across cells.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
sc.pp.highly_variable_genes(adata, n_top_genes = 2000)
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata.var
```

```python slideshow={"slide_type": "subslide"} tags=[]
sc.pl.highly_variable_genes(adata)
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Filter highly variable genes
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
We finally filter for the highly variable genes noting the change in `n_vars`.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
adata
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata = adata[:, adata.var.highly_variable]
```

```python slideshow={"slide_type": "fragment"} tags=[]
adata
```

<!-- #region {"slideshow": {"slide_type": "slide"}, "tags": []} -->
### Principle components analysis
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
#### Regression and rescaling
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
We would like to avoid interpreting differences that derive from differences in total counts, percentage of mitochondrial counts, and percentage of ribosomal counts. We can reduce the impact of these differences using the [`scanpy` regress out function](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.regress_out.html).

<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt', 'pct_counts_ribo'])
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
We also want to normalize the data [scaling](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.scale.html) it to zero mean and unit variance while clipping any remaining outliers above 10.
<!-- #endregion -->

```python slideshow={"slide_type": "fragment"} tags=[]
sc.pp.scale(adata, max_value=10)
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
#### Perform PCA
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
PCA can be performed with the [`scanpy` PCA tool](https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.pca.html).
<!-- #endregion -->

```python
sc.tl.pca(adata, svd_solver='arpack')
```

```python
sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50)
```

### Clustering


The UMAP algorithm contains a method for estimating distances between cells to which an interface is provided by [scanpy neighbors](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.neighbors.html#scanpy.pp.neighbors).

```python
sc.pp.neighbors(adata, n_pcs = 30)
```

This function adds `distances` and `connectivities` as pairwise relationships among the cells in `adata.obsp`.

```python
adata
```

```python
adata.obsp
```

We can now plot the UMAP embedding, but note there are no defined clusters of cells even though they may be hypothesized visually.

```python
sc.tl.umap(adata)
```

```python
sc.pl.umap(adata)
```

The [`scanpy` interface to the Leiden algorithm](https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.leiden.html) can be used to cluster cells into subgroups.

```python
sc.tl.leiden(adata, resolution = 0.5)
```

```python
adata.obs
```

We can now color the UMAP embedding by the Louvain clusters.

```python
sc.pl.umap(adata, color=['leiden'])
```

<!-- #region {"slideshow": {"slide_type": "subslide"}, "tags": []} -->
### Save/load checkpoint 
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Save current data with cluster information to file.
<!-- #endregion -->

```python tags=[]
adata.write(f"{output_directory}/adata_clustered.h5ad", compression="gzip")
```

```python
adata
```

<!-- #region {"slideshow": {"slide_type": "fragment"}, "tags": []} -->
Variables can be reloaded if necessary.
<!-- #endregion -->

```python
bdata_clustered = None
bdata_clustered = sc.read_h5ad(f"{output_directory}/adata_clustered.h5ad")
bdata_clustered
```

## Sample integration


### Batch processing

<!-- #region -->
See the command line program [scpreprocess](./sc/scpreprocess.py). It can be installed with

```sh
$ pip install -e .
```

and run, for example

```sh
$ scpreprocess  --input-directory="data/GSE171524/supplementary/"
```

Where the contents of input-directory was, in this case, created by running [download_geo_data.sh](./data/download_geo_data.sh) on `GSE171524`. In this case, `scpreprocess` requires about 30GB of RAM and about 30 - 60 minutes to complete with 1 NVIDIA Tesla T4 GPU.
<!-- #endregion -->

### Filter genes

```python tags=[]
adata = None
adata = sc.read_h5ad(f"{output_directory}/adata_combined.h5ad")
adata
```

```python tags=[]
sc.pp.filter_genes(adata, min_cells = 10)
```

```python
adata.X = csr_matrix(adata.X)
```

```python
adata.write(f"{output_directory}/adata_combined_filtered.h5ad", compression="gzip")
```

```python
adata = None
adata = sc.read_h5ad(f"{output_directory}/adata_combined_filtered.h5ad")
adata
```

```python
adata.obs.groupby('Sample').count()
```

```python
sc.pp.filter_genes(adata, min_cells = 100)
```

```python
adata
```

```python
adata.layers['counts'] = adata.X.copy()
```

```python
adata
```

### Normalize

```python
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
adata.raw = adata
```

```python
adata.obs.head()
```

```python
scvi.model.SCVI.setup_anndata(adata, layer = "counts",
                             categorical_covariate_keys=["Sample"], # batch, technology
                             continuous_covariate_keys=['pct_counts_mt', 'total_counts', 'pct_counts_ribo'])
```

```python
model = scvi.model.SCVI(adata)
```

```python
model.train()
```

```python tags=[]
# pickle.dump([model], open(f"{output_directory}/model_combined.p", "wb"))
model.save(f"{output_directory}/model_combined/")
# model = scvi.model.SCVI.load(f"{output_directory}/model_combined/")
```

```python
adata.obsm['X_scVI'] = model.get_latent_representation()
```

```python
adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)
```

### Cluster

```python
sc.pp.neighbors(adata, use_rep = 'X_scVI')
```

```python tags=[]
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 0.5)
```

```python tags=[]
sc.pl.umap(adata, color = ['leiden'], frameon = True)
```

```python
sc.pl.umap(adata, color = ['Sample'], frameon = True)
```

### Save/load checkpoint

```python tags=[]
adata.write(f"{output_directory}/adata_integrated.h5ad", compression="gzip")
```

```python
adata = None
adata = sc.read_h5ad(f"{output_directory}/adata_integrated.h5ad")
adata
```

## Marker identification


Increase cluster resolution from 0.5 to [default value of 1](https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.leiden.html#scanpy.tl.leiden).

```python
sc.tl.leiden(adata, resolution = 1)
```

```python jupyter={"outputs_hidden": true} tags=[]
sc.tl.rank_genes_groups(adata, 'leiden')
```

Plot genes that distinguish each cluster.

```python
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
```

Filter markers by adjusted p-value and log fold change.

```python
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers
```




Estimate differential expression for leiden clusters.

```python
markers_scvi = model.differential_expression(groupby = 'leiden')
markers_scvi
```

Filter differentially expressed genes by FDR and log fold change.

```python
markers_scvi = markers_scvi[(markers_scvi['is_de_fdr_0.05']) & (markers_scvi.lfc_mean > .5)]
markers_scvi
```

```python tags=[]
sc.pl.umap(adata, color = ['leiden'], frameon = False, legend_loc = "on data")
```

```python
sc.pl.umap(adata, color = ['PTPRC', 'CD3E', 'CD4'], frameon = False, layer = 'scvi_normalized', vmax = 5)
sc.pl.umap(adata, color = ['PTPRC', 'CD3E', 'CD8A'], frameon = False, layer = 'scvi_normalized', vmax = 5)

sc.pl.umap(adata, color = ['EPCAM', 'MUC1'], frameon = False, layer = 'scvi_normalized', vmax = 5)
sc.pl.umap(adata, color = ['EPCAM', 'MUC1'], frameon = False, layer = 'scvi_normalized', vmax = 5)
```
