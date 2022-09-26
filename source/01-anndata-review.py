# ---
# jupyter:
#   anaconda-cloud: {}
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

# %% [markdown]
# # AnnData review

# %% [markdown]
# This is a modified version of the introductory tutorial distributed with the AnnData package.
#
# **Authors:** [Adam Gayoso](https://twitter.com/adamgayoso), [Alex Wolf](https://twitter.com/falexwolf)
#
# This tutorial is based on blog posts by [Adam in 2021](https://adamgayoso.com/posts/ten_min_to_adata/) and [Alex in 2017](https://falexwolf.me/2017/introducing-anndata/).
#
# This tutorial introduces basic properties of the central object, [AnnData](http://anndata.readthedocs.io/en/latest/anndata.AnnData.html) ("Annotated Data").
#
# `AnnData` is specifically designed for matrix-like data. By this we mean that we have $n$ observations, each of which can be represented as $d$-dimensional vectors, where each dimension corresponds to a variable or feature. Both the rows and columns of this $n \times d$ matrix are special in the sense that they are indexed.
#
# For instance, in scRNA-seq data, each row corresponds to a cell with a barcode, and each column corresponds to a gene with a gene id. Furthermore, for each cell and each gene we might have additional metadata, like (1) donor information for each cell, or (2) alternative gene symbols for each gene. Finally, we might have other unstructured metadata like color palletes to use for plotting. Without going into every fancy Python-based data structure, we think that still today no other alternative really exists that:
#
# * Handles sparsity
# * Handles unstructured data
# * Handles observation- and feature-level metadata
# * Is user-friendly

# %% {"jupyter": {"outputs_hidden": false}, "tags": []}
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
print(ad.__version__)

# %% [markdown]
# ### Utility functions

# %% {"tags": []}
from pprint import pprint
from inspect import getmembers
from types import FunctionType


# %% {"tags": []}
def attributes(obj):
    disallowed_names = {
        name for name, value in getmembers(type(obj)) if isinstance(value, FunctionType)
    }
    return {
        name: getattr(obj, name)
        for name in dir(obj)
        if name[0] != "_" and name not in disallowed_names and hasattr(obj, name)
    }


def print_attributes(obj):
    pprint(attributes(obj))


# %% [markdown]
# ## Initializing AnnData

# %% [markdown]
# Let's start by building a basic AnnData object with some sparse count information, perhaps representing gene expression counts.

# %% {"tags": []}
np.random.poisson(0.5, size=(4, 9))

# %% {"tags": []}
counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
# counts = csr_matrix(np.random.poisson(1, size=(4, 9)), dtype=np.float32)

# %% {"tags": []}
counts

# %% {"tags": []}
adata = ad.AnnData(counts)
adata

# %% {"tags": []}
print_attributes(adata)

# %% [markdown]
# We see that AnnData provides a representation with summary stastics of the data The initial data we passed are accessible as a sparse matrix using `adata.X`.

# %% {"tags": []}
adata.X

# %% [markdown]
# Now, we provide the index to both the `obs` and `var` axes using `.obs_names` (resp. `.var_names`).

# %% {"tags": []}
print(adata.n_obs)
print(adata.n_vars)

# %% {"tags": []}
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]

# %% {"tags": []}
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]

# %% {"tags": []}
print(adata.obs_names[:10])

# %% {"tags": []}
print(adata.var_names[:10])

# %% {"tags": []}
print_attributes(adata)

# %% [markdown]
# ### Subsetting AnnData

# %% [markdown]
# These index values can be used to subset the AnnData, which provides a view of the AnnData object. We can imagine this to be useful to subset the AnnData to particular cell types or gene modules of interest. The rules for subsetting AnnData are quite similar to that of a Pandas DataFrame. You can use values in the `obs/var_names`, boolean masks, or cell index integers.

# %%
adata[["Cell_1", "Cell_10"], ["Gene_5", "Gene_1900"]]

# %% [markdown]
# ## Adding aligned metadata
#
# ### Observation/Variable level

# %% [markdown]
# So we have the core of our object and now we'd like to add metadata at both the observation and variable levels. This is pretty simple with AnnData, both `adata.obs` and `adata.var` are Pandas DataFrames.

# %%
ct = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
adata.obs["cell_type"] = pd.Categorical(ct)  # Categoricals are preferred for efficiency
adata.obs

# %% [markdown]
# We can also see now that the AnnData representation has been updated:

# %%
adata

# %% [markdown]
# ### Subsetting using metadata

# %% [markdown]
# We can also subset the AnnData using these randomly generated cell types:

# %%
bdata = adata[adata.obs.cell_type == "B"]
bdata

# %% [markdown]
# ## Observation/variable-level matrices

# %% [markdown]
# We might also have metadata at either level that has many dimensions to it, such as a UMAP embedding of the data. For this type of metadata, AnnData has the `.obsm/.varm` attributes. We use keys to identify the different matrices we insert. The restriction of `.obsm/.varm` are that `.obsm` matrices must length equal to the number of observations as `.n_obs` and `.varm` matrices must length equal to `.n_vars`. They can each independently have different number of dimensions.
#
# Let's start with a randomly generated matrix that we can interpret as a UMAP embedding of the data we'd like to store, as well as some random gene-level metadata:

# %%
adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))
adata.varm["gene_stuff"] = np.random.normal(0, 1, size=(adata.n_vars, 5))
adata.obsm

# %% [markdown]
# Again, the AnnData representation is updated.

# %%
adata

# %% [markdown]
# A few more notes about `.obsm/.varm`
#
# 1. The "array-like" metadata can originate from a Pandas DataFrame, scipy sparse matrix, or numpy dense array.
# 2. When using scanpy, their values (columns) are not easily plotted, where instead items from `.obs` are easily plotted on, e.g., UMAP plots.

# %% [markdown]
# ## Unstructured metadata

# %% [markdown]
# AnnData has `.uns`, which allows for any unstructured metadata. This can be anything, like a list or a dictionary with some general information that was useful in the analysis of our data.

# %%
adata.uns["random"] = [1, 2, 3]
adata.uns

# %% [markdown]
# ## Layers

# %% [markdown]
# Finally, we may have different forms of our original core data, perhaps one that is normalized and one that is not. These can be stored in different layers in AnnData. For example, let's log transform the original data and store it in a layer:

# %%
adata.layers["log_transformed"] = np.log1p(adata.X)
adata

# %% [markdown]
# ## Conversion to DataFrames

# %% [markdown]
# We can also ask AnnData to return us a DataFrame from one of the layers:

# %%
adata.to_df(layer="log_transformed")

# %% [markdown]
# We see that the `.obs_names/.var_names` are used in the creation of this Pandas object.

# %% [markdown]
# ## Writing the results to disk

# %% [markdown]
# `AnnData` comes with its own persistent HDF5-based file format: `h5ad`. If string columns with small number of categories aren't yet categoricals, `AnnData` will auto-transform to categoricals.

# %%
adata.write('my_results.h5ad', compression="gzip")

# %% {"jupyter": {"outputs_hidden": false}}
# !h5ls 'my_results.h5ad'

# %% [markdown]
# ## Wrapping up the introduction

# %% [markdown]
# AnnData has become the standard for single-cell analysis in Python and for good reason -- it's straightforward to use and faciliatates more reproducible analyses with it's key-based storage. It's even becoming easier to convert to the popular R-based formats for single-cell analysis.
#
# Keep reading on to better understand "views", on-disk backing, and other details.

# %% [markdown]
# ## Views and copies

# %% [markdown]
# For the fun of it, let's look at another metadata use case. Imagine that the observations come from instruments characterizing 10 readouts in a multi-year study with samples taken from different subjects at different sites. We'd typically get that information in some format and then store it in a DataFrame:

# %% {"jupyter": {"outputs_hidden": false}}
obs_meta = pd.DataFrame({
        'time_yr': np.random.choice([0, 2, 4, 8], adata.n_obs),
        'subject_id': np.random.choice(['subject 1', 'subject 2', 'subject 4', 'subject 8'], adata.n_obs),
        'instrument_type': np.random.choice(['type a', 'type b'], adata.n_obs),
        'site': np.random.choice(['site x', 'site y'], adata.n_obs),
    },
    index=adata.obs.index,    # these are the same IDs of observations as above!
)

# %% [markdown]
# This is how we join the readout data with the metadata. Of course, the first argument of the following call for `X` could also just be a DataFrame.

# %% {"jupyter": {"outputs_hidden": false}}
adata = ad.AnnData(adata.X, obs=obs_meta, var=adata.var)

# %% [markdown]
# Now we again have a single data container that keeps track of everything.

# %% {"jupyter": {"outputs_hidden": false}}
print(adata)

# %% [markdown]
# Subsetting the joint data matrix can be important to focus on subsets of variables or observations, or to define train-test splits for a machine learning model.

# %% [markdown]
# <div class="alert alert-info">
#
# **Note**
#
# Similar to numpy arrays, AnnData objects can either hold actual data or reference another `AnnData` object. In the later case, they are referred to as "view".
#
# Subsetting AnnData objects always returns views, which has two advantages:
#
# - no new memory is allocated
# - it is possible to modify the underlying AnnData object
#
# You can get an actual AnnData object from a view by calling `.copy()` on the view. Usually, this is not necessary, as any modification of elements of a view (calling `.[]` on an attribute of the view) internally calls `.copy()` and makes the view an AnnData object that holds actual data. See the example below.
#     
# </div>

# %% {"jupyter": {"outputs_hidden": false}}
adata

# %% [markdown]
# Get access to the first 5 rows for two variables.

# %% [markdown]
# <div class="alert alert-info">
#
# **Note**
#
# Indexing into AnnData will assume that integer arguments to `[]` behave like `.iloc` in pandas, whereas string arguments behave like `.loc`. `AnnData` always assumes string indices.
#     
# </div> 

# %% {"jupyter": {"outputs_hidden": false}}
adata[:5, ['Gene_1', 'Gene_3']]

# %% [markdown]
# This is a view! If we want an `AnnData` that holds the data in memory, let's call `.copy()`

# %%
adata_subset = adata[:5, ['Gene_1', 'Gene_3']].copy()

# %% [markdown]
# For a view, we can also set the first 3 elements of a column.

# %% {"jupyter": {"outputs_hidden": false}}
print(adata[:3, 'Gene_1'].X.toarray().tolist())
adata[:3, 'Gene_1'].X = [0, 0, 0]
print(adata[:3, 'Gene_1'].X.toarray().tolist())

# %% [markdown]
# If you try to access parts of a view of an AnnData, the content will be auto-copied and a data-storing object will be generated.

# %% {"jupyter": {"outputs_hidden": false}}
adata_subset = adata[:3, ['Gene_1', 'Gene_2']]

# %% {"jupyter": {"outputs_hidden": false}}
adata_subset

# %% {"jupyter": {"outputs_hidden": false}}
adata_subset.obs['foo'] = range(3)

# %% [markdown]
# Now `adata_subset` stores the actual data and is no longer just a reference to `adata`.

# %% {"jupyter": {"outputs_hidden": false}}
adata_subset

# %% [markdown]
# Evidently, you can use all of pandas to slice with sequences or boolean indices.

# %% {"jupyter": {"outputs_hidden": false}}
adata[adata.obs.time_yr.isin([2, 4])].obs.head()

# %% [markdown]
# ## Partial reading of large data

# %% [markdown]
# If a single `.h5ad` is very large, you can partially read it into memory by using backed mode:

# %%
adata = ad.read('my_results.h5ad', backed='r')

# %%
adata.isbacked

# %% [markdown]
# If you do this, you'll need to remember that the `AnnData` object has an open connection to the file used for reading:

# %%
adata.filename

# %% [markdown]
# As we're using it in read-only mode, we can't damage anything. To proceed with this tutorial, we still need to explicitly close it: 

# %%
adata.file.close()

# %% [markdown]
# As usual, you should rather use `with` statements to avoid dangling open files (up-coming feature).

# %% [markdown]
# Manipulating the object on disk is possible, but experimental for sparse data. Hence, we leave it out of this tutorial.
