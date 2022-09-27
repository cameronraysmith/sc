# ---
# jupyter:
#   celltoolbar: Slideshow
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,md,py:percent
#     notebook_metadata_filter: all
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
#   language_info:
#     codemirror_mode:
#       name: ipython
#       version: 3
#     file_extension: .py
#     mimetype: text/x-python
#     name: python
#     nbconvert_exporter: python
#     pygments_lexer: ipython3
#     version: 3.10.7
#   rise:
#     scroll: true
#     theme: black
#   toc-autonumbering: true
#   toc-showcode: false
#   toc-showmarkdowntxt: false
#   widgets:
#     application/vnd.jupyter.widget-state+json:
#       state: {}
#       version_major: 2
#       version_minor: 0
# ---

# %% [markdown] {"slideshow": {"slide_type": "slide"}, "tags": []}
# # Basic single-cell analysis

# %% [markdown] {"slideshow": {"slide_type": "subslide"}, "tags": []}
# ## Overview

# %% [markdown] {"slideshow": {"slide_type": "fragment"}, "tags": []}
# This notebook follows the tutorial by [mousepixels/sanbomics](https://github.com/mousepixels/sanbomics/blob/main/single_cell_analysis_complete_class.ipynb), which has an accompanying [screencast](https://youtu.be/uvyG9yLuNSE?t=319).
#
# Analysis is illustrated with single-nucleus RNA sequencing data from the following paper <cite data-cite="Melms2021-bj">Melms et al. (2021)</cite>
#
# > Melms JC, Biermann J, Huang H, Wang Y, Nair A, Tagore S, et al.
# A molecular single-cell lung atlas of lethal COVID-19.
# Nature. 2021;595: 114–119. [doi:10.1038/s41586-021-03569-0](https://doi.org/10.1038/s41586-021-03569-0)
#
# This paper examined 116,000 nuclei from the lungs of nineteen patients who underwent autopsy following death in association with COVID-19. Findings reported in the abstract of the paper include:
#
# 1. activated monocyte-derived macrophages and alveolar macrophages
# 1. impaired T cell activation
# 1. monocyte/macrophage-derived interleukin-1β and epithelial cell-derived interleukin-6
# 1. alveolar type 2 cells adopted an inflammation-associated transient progenitor cell state and failed to undergo full transition into alveolar type 1 cells
# 1. expansion of CTHRC1+ pathological fibroblasts
# 1. protein activity and ligand–receptor interactions suggest putative drug targets
#
# This notebook makes extensive use of <cite data-cite="Wolf2018-nu">Wolf et al. (2018)</cite> and <cite data-cite="Lopez2018-em">Lopez et al. (2018)</cite> including updates that have been made to the underlying software packages, [scanpy](https://github.com/scverse/scanpy) and [scvi-tools](https://github.com/scverse/scvi-tools), since their initial publication.

# %% [markdown] {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"slide\"} jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "slide"}, "tags": []}
# ## Setup

# %% [markdown] {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"subslide\"} jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "subslide"}, "tags": []}
# ### Import libraries

# %% {"tags": [], "slideshow": {"slide_type": "fragment"}}
from inspect import getmembers
from pprint import pprint
from types import FunctionType

import scanpy as sc

# %% [markdown] {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"subslide\"} jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "subslide"}, "tags": []}
# ### Setup plotting

# %% {"slideshow": {"slide_type": "fragment"}, "tags": []}
import matplotlib.font_manager
import matplotlib.pyplot as plt

# import matplotlib_inline

# %% {"slideshow": {"slide_type": "fragment"}, "tags": []}
# fonts_path = "/usr/share/texmf/fonts/opentype/public/lm/" #ubuntu
# fonts_path = "~/Library/Fonts/" # macos
fonts_path = "/usr/share/fonts/OTF/"  # arch
# user_path = "$HOME/" # user
# fonts_path = user_path + "fonts/latinmodern/opentype/public/lm/"  # home
matplotlib.font_manager.fontManager.addfont(fonts_path + "lmsans10-regular.otf")
matplotlib.font_manager.fontManager.addfont(fonts_path + "lmroman10-regular.otf")

# %% {"slideshow": {"slide_type": "fragment"}, "tags": []}
# https://stackoverflow.com/a/36622238/446907
# %config InlineBackend.figure_formats = ['svg']

# %% {"slideshow": {"slide_type": "fragment"}, "tags": []}
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


# %% [markdown] {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"subslide\"} jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "subslide"}, "tags": []}
# ### Utility functions

# %% {"tags": [], "slideshow": {"slide_type": "fragment"}}
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


# %% [markdown] {"slideshow": {"slide_type": "slide"}, "tags": []}
# ## Import data

# %% [markdown] {"slideshow": {"slide_type": "fragment"}, "tags": []}
# Here we review how the data were downloaded, and proceed to import and inspect the data.

# %% [markdown] {"slideshow": {"slide_type": "subslide"}, "tags": []}
# ### Data download

# %% [markdown] {"slideshow": {"slide_type": "fragment"}, "tags": []}
# Data with GEO accession [GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524) was downloaded using [./data/download_geo_data.sh](./data/download_geo_data.sh) with parameters
#
# ```bash
# ./download_geo_data.sh \
#        -a GSE132771 \
#        -f 'ftp.*RAW.*' \
#        -j '..|.supplementary_files?|..|.url?|select(length>0)'
# ```
#
# A skeleton of this script that may work in this case is
#
# ```bash
# !/usr/bin/env bash
#
# #-- debugging (comment to reduce stderr output)
# #-- https://wiki.bash-hackers.org/scripting/debuggingtips
# export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
# set -o xtrace
#
# # get metadata
# # Melms JC, Biermann J, Huang H, Wang Y, Nair A, Tagore S, et al.
# # A molecular single-cell lung atlas of lethal COVID-19.
# # Nature. 2021;595: 114–119. doi:10.1038/s41586-021-03569-0
# # GSE171524
# ffq -l 1 -o GSE171524.json GSE171524
#
# # download raw data
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171524/suppl/GSE171524_RAW.tar
#
# # list contents
# tar -tvf GSE171524_RAW.tar
#
# # untar
# mkdir -p GSE171524 && \
# tar -xvf GSE171524_RAW.tar -C GSE171524
# ```

# %% [markdown] {"slideshow": {"slide_type": "subslide"}, "tags": []}
# ### Data load

# %% {"slideshow": {"slide_type": "fragment"}, "tags": []}
adata = None
adata = sc.read_csv("data/GSE171524/supplementary/GSM5226574_C51ctr_raw_counts.csv.gz").T
adata

# %% [markdown] {"slideshow": {"slide_type": "fragment"}, "tags": []}
# Note the `scanpy.read_csv` function accepts gzipped files.

# %% [markdown] {"slideshow": {"slide_type": "subslide"}, "tags": []}
# ### Data properties

# %% {"slideshow": {"slide_type": "fragment"}, "tags": []}
type(adata)

# %% {"tags": []}
type(adata.T)

# %% {"tags": [], "slideshow": {"slide_type": "fragment"}}
print_attributes(adata)

# %% {"tags": []}
adata.obs

# %% [markdown]
# Gene names are saved 

# %% {"tags": []}
adata.var

# %% {"tags": []}
adata.obs_names

# %% {"tags": []}
adata.var_names

# %% [markdown]
# There are two layers corresponding to spliced and unspliced transcripts respectively.

# %% {"tags": []}
adata.layers['spliced']

# %% {"tags": []}
adata.layers['unspliced']

# %% [markdown]
# PCA and UMAP have retained 50 and 2 dimensions respectively.

# %% {"tags": []}
print(adata.obsm)
print(adata.obsm['X_pca'].shape)
print(adata.obsm)
print(adata.obsm['X_umap'].shape)

# %% {"tags": []}
print(adata.varm)

# %% {"tags": []}
print(adata.obsp)
print(adata.obsp['distances'].shape)
print(adata.obsp)
print(adata.obsp['connectivities'].shape)

# %% {"tags": []}
print(adata.varp)

# %% [markdown] {"slideshow": {"slide_type": "fragment"}, "tags": []}
# The data appears to contain reads mapped to 34546 RNA molecule-associated features and 6099 cell-associated barcodes.

# %%
