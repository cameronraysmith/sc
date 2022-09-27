# ---
# jupyter:
#   celltoolbar: Slideshow
#   hide_code_all_hidden: false
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
#   nbpresent:
#     slides: {}
#     themes:
#       default: 15f06427-3b2a-49bc-bf86-daf4b06df6ed
#       theme:
#         15f06427-3b2a-49bc-bf86-daf4b06df6ed:
#           backgrounds:
#             dc7afa04-bf90-40b1-82a5-726e3cff5267:
#               background-color: 31af15d2-7e15-44c5-ab5e-e04b16a89eff
#               id: dc7afa04-bf90-40b1-82a5-726e3cff5267
#           id: 15f06427-3b2a-49bc-bf86-daf4b06df6ed
#           palette:
#             19cc588f-0593-49c9-9f4b-e4d7cc113b1c:
#               id: 19cc588f-0593-49c9-9f4b-e4d7cc113b1c
#               rgb:
#               - 252
#               - 252
#               - 252
#             31af15d2-7e15-44c5-ab5e-e04b16a89eff:
#               id: 31af15d2-7e15-44c5-ab5e-e04b16a89eff
#               rgb:
#               - 68
#               - 68
#               - 68
#             50f92c45-a630-455b-aec3-788680ec7410:
#               id: 50f92c45-a630-455b-aec3-788680ec7410
#               rgb:
#               - 197
#               - 226
#               - 245
#             c5cc3653-2ee1-402a-aba2-7caae1da4f6c:
#               id: c5cc3653-2ee1-402a-aba2-7caae1da4f6c
#               rgb:
#               - 43
#               - 126
#               - 184
#             efa7f048-9acb-414c-8b04-a26811511a21:
#               id: efa7f048-9acb-414c-8b04-a26811511a21
#               rgb:
#               - 25.118061674008803
#               - 73.60176211453744
#               - 107.4819383259912
#           rules:
#             a:
#               color: 19cc588f-0593-49c9-9f4b-e4d7cc113b1c
#             blockquote:
#               color: 50f92c45-a630-455b-aec3-788680ec7410
#               font-size: 3
#             code:
#               font-family: Anonymous Pro
#             h1:
#               color: 19cc588f-0593-49c9-9f4b-e4d7cc113b1c
#               font-family: Merriweather
#               font-size: 8
#             h2:
#               color: 19cc588f-0593-49c9-9f4b-e4d7cc113b1c
#               font-family: Merriweather
#               font-size: 6
#             h3:
#               color: 50f92c45-a630-455b-aec3-788680ec7410
#               font-family: Lato
#               font-size: 5.5
#             h4:
#               color: c5cc3653-2ee1-402a-aba2-7caae1da4f6c
#               font-family: Lato
#               font-size: 5
#             h5:
#               font-family: Lato
#             h6:
#               font-family: Lato
#             h7:
#               font-family: Lato
#             li:
#               color: 50f92c45-a630-455b-aec3-788680ec7410
#               font-size: 3.25
#             pre:
#               font-family: Anonymous Pro
#               font-size: 4
#           text-base:
#             color: 19cc588f-0593-49c9-9f4b-e4d7cc113b1c
#             font-family: Lato
#             font-size: 4
#   rise:
#     scroll: true
#     theme: black
#   toc-autonumbering: true
#   toc-showcode: false
#   toc-showmarkdowntxt: false
#   varInspector:
#     cols:
#       lenName: 16
#       lenType: 16
#       lenVar: 40
#     kernels_config:
#       python:
#         delete_cmd_postfix: ''
#         delete_cmd_prefix: 'del '
#         library: var_list.py
#         varRefreshCmd: print(var_dic_list())
#       r:
#         delete_cmd_postfix: ') '
#         delete_cmd_prefix: rm(
#         library: var_list.r
#         varRefreshCmd: 'cat(var_dic_list()) '
#     types_to_exclude:
#     - module
#     - function
#     - builtin_function_or_method
#     - instance
#     - _Feature
#     window_display: false
#   widgets:
#     application/vnd.jupyter.widget-state+json:
#       state: {}
#       version_major: 2
#       version_minor: 0
# ---

# %% [markdown] {"hideCode": false, "hidePrompt": false, "slideshow": {"slide_type": "slide"}, "tags": []}
# # RNA Velocity review

# %% [markdown] {"slideshow": {"slide_type": "fragment"}, "tags": []}
# Here we review basic RNA velocity analysis.
#
# This is a modified version of the [scvelo basic velocity analysis](https://nbviewer.jupyter.org/github/theislab/scvelo_notebooks/blob/master/VelocityBasics.ipynb) tutorial.
#
# The primary data set is development of the endocrine pancreas, which demonstrates lineage commitment to four distinct fates: α, β, δ and ε-cells. <br/> 
# See [here](https://scvelo.readthedocs.io/scvelo.datasets.pancreas.html) for more details.
#
#

# %% [markdown] {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"slide\"} jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "slide"}, "tags": []}
# ## Setup

# %% [markdown] {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"subslide\"} jp-MarkdownHeadingCollapsed=true slideshow={\"slide_type\": \"subslide\"} tags=[] jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "subslide"}, "tags": []}
# ### Import libraries

# %% {"slideshow": {"slide_type": "fragment"}, "tags": []}
# update to the latest version, if not done yet.
# # !pip install scvelo --upgrade --quiet

# %% {"tags": [], "slideshow": {"slide_type": "fragment"}}
from inspect import getmembers
from pprint import pprint
from types import FunctionType

import scanpy as sc
import scvelo as scv

# %% {"hideCode": false, "hidePrompt": false, "slideshow": {"slide_type": "fragment"}, "tags": []}
scv.logging.print_version()

# %% {"hideCode": false, "hidePrompt": false, "slideshow": {"slide_type": "fragment"}, "tags": []}
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


# %% [markdown] {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"subslide\"} jp-MarkdownHeadingCollapsed=true slideshow={\"slide_type\": \"subslide\"} tags=[] jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "subslide"}, "tags": []}
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


# %% [markdown] {"incorrectly_encoded_metadata": "tags=[] slideshow={\"slide_type\": \"subslide\"} jp-MarkdownHeadingCollapsed=true slideshow={\"slide_type\": \"subslide\"} tags=[] jp-MarkdownHeadingCollapsed=true", "slideshow": {"slide_type": "subslide"}, "tags": []}
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


# %% [markdown] {"hideCode": false, "hidePrompt": false, "slideshow": {"slide_type": "slide"}, "tags": []}
# ## Load data

# %% [markdown] {"slideshow": {"slide_type": "fragment"}, "tags": []}
# * load an example dataset into an `AnnData` object
#     - the data set used here contains precomputed spliced/unspliced reads
#     - spliced/unspliced reads can be computed with `velocyto.py` and `kb-python`
# * view the detailed structure of the object
# * estimate relative proportions of spliced/unspliced reads in the data

# %% [markdown] {"slideshow": {"slide_type": "subslide"}, "tags": []}
# ### Load AnnData object

# %% [markdown] {"hideCode": false, "hidePrompt": false, "slideshow": {"slide_type": "fragment"}, "tags": []}
# The analysis is based on the in-built [pancreas data](https://scvelo.readthedocs.io/scvelo.datasets.pancreas).<br/>
# To run velocity analysis on your own data, read your file (loom, h5ad, csv …) to an AnnData object with `adata = scv.read('path/file.loom', cache=True)`. If you want to merge your loom file into an already existing AnnData object, use `scv.utils.merge(adata, adata_loom)`.

# %% {"hideCode": false, "hidePrompt": false, "slideshow": {"slide_type": "fragment"}, "tags": []}
adata = scv.datasets.pancreas()
adata

# %% [markdown] {"hideCode": false, "hidePrompt": false, "slideshow": {"slide_type": "fragment"}, "tags": []}
# scVelo is based on `adata`, an object that stores a data matrix `adata.X`, annotation of observations `adata.obs`, variables `adata.var`, and unstructured annotations `adata.uns`. Names of observations and variables can be accessed via `adata.obs_names` and `adata.var_names`, respectively. AnnData objects can be sliced like dataframes, for example, `adata_subset = adata[:, list_of_gene_names]`. For more details, see the [anndata docs](https://anndata.readthedocs.io).
#

# %% [markdown] {"slideshow": {"slide_type": "subslide"}, "tags": []}
# ### Detailed view of adata

# %% [markdown]
# Observations are indexed by barcodes that are intended to correspond to cells after preprocessing. There are four components to the annotation vector. Two of these are fine- and coarse-grained clusters into cell types. The others are scores `S_score` and `G2M_score` that assess expression of G2/M and S phase marker genes.

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

# %% {"slideshow": {"slide_type": "fragment"}, "tags": []}
print_attributes(adata)

# %% [markdown] {"slideshow": {"slide_type": "subslide"}, "tags": []}
# ### Show spliced/unspliced

# %% [markdown] {"hideCode": false, "hidePrompt": false, "slideshow": {"slide_type": "fragment"}, "tags": []}
# Here, the proportions of spliced/unspliced counts are displayed. Depending on the protocol used (Drop-Seq, Smart-Seq), we typically have between 10%-25% of unspliced molecules containing intronic sequences. We also advice you to examine the variations on cluster level to verify consistency in splicing efficiency. Here, we find variations as expected, with slightly lower unspliced proportions at cycling ductal cells, then higher proportion at cell fate commitment in Ngn3-high and Pre-endocrine cells where many genes start to be transcribed.

# %% {"hideCode": false, "hidePrompt": false, "slideshow": {"slide_type": "fragment"}, "tags": []}
scv.pl.proportions(adata)

# %% [markdown] {"hideCode": false, "hidePrompt": false, "slideshow": {"slide_type": "slide"}, "tags": []}
# ## Preprocess the Data

# %% [markdown] {"hideCode": false, "hidePrompt": false, "slideshow": {"slide_type": "fragment"}, "tags": []}
# Preprocessing requisites consist of **gene selection** by detection (with a minimum number of counts) and high variability (dispersion), **normalizing** every cell by its total size and **logarithmizing** X. Filtering and normalization is applied in the same vein to spliced/unspliced counts and X. Logarithmizing is only applied to X. If X is already preprocessed from former analysis, it will not be touched.
#
# All of this is summarized in a single function `scv.pp.filter_and_normalize`, which essentially runs the following:
#
# ```
# scv.pp.filter_genes(adata, min_shared_counts=20)
# scv.pp.normalize_per_cell(adata)
# scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
# scv.pp.log1p(adata)
# ```
#
# Further, we need the first and second order moments (means and uncentered variances) computed among nearest neighbors in PCA space, summarized in `scv.pp.moments`, which internally computes `scv.pp.pca` and `scv.pp.neighbors`.
# First order is needed for deterministic velocity estimation, while stochastic estimation also requires second order moments.

# %% {"hideCode": false, "hidePrompt": false}
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


# %% [markdown] {"hideCode": false, "hidePrompt": false}
# Further preprocessing (such as batch effect correction) may be used to remove unwanted sources of variability. See the [best practices](https://doi.org/10.15252/msb.20188746) for further details. Note, that any additional preprocessing step only affects X and is not applied to spliced/unspliced counts.

# %% [markdown]
# ## Compute RNA velocity

# %% [markdown] {"hideCode": false, "hidePrompt": false}
# ### Estimate RNA velocity

# %% [markdown] {"hideCode": false, "hidePrompt": false}
# Velocities are vectors in gene expression space and represent the direction and speed of movement of the individual cells. The velocities are obtained by modeling transcriptional dynamics of splicing kinetics, either stochastically (default) or deterministically (by setting ``mode='deterministic'``). For each gene, a steady-state-ratio of pre-mature (unspliced) and mature (spliced) mRNA counts is fitted, which constitutes a constant transcriptional state. Velocities are then obtained as residuals from this ratio. Positive velocity indicates that a gene is up-regulated, which occurs for cells that show higher abundance of unspliced mRNA for that gene than expected in steady state. Conversely, negative velocity indicates that a gene is down-regulated.
#
# The solution to the full dynamical model is obtained by setting ``mode='dynamical'``, which requires to run
# ``scv.tl.recover_dynamics(adata)`` beforehand. We will elaborate more on the dynamical model in the next tutorial.

# %% {"hideCode": false, "hidePrompt": false}
scv.tl.velocity(adata)

# %% [markdown] {"hideCode": false, "hidePrompt": false}
# The computed velocities are stored in ``adata.layers`` just like the count matrices. 
#
# The combination of velocities across genes can then be used to estimate the future state of an individual cell. In order to project the velocities into a lower-dimensional embedding, transition probabilities of cell-to-cell transitions are estimated. That is, for each velocity vector we find the likely cell transitions that are accordance with that direction. The transition probabilities are computed using cosine correlation between the potential cell-to-cell transitions and the velocity vector, and are stored in a matrix denoted as velocity graph. The resulting velocity graph has dimension $n_{obs} \times n_{obs}$ and summarizes the possible cell state changes that are well explained through the velocity vectors (for runtime speedup it can also be computed on reduced PCA space by setting `approx=True`).

# %% {"hideCode": false, "hidePrompt": false}
scv.tl.velocity_graph(adata)

# %% [markdown] {"hideCode": false, "hidePrompt": false}
# For a variety of applications, the velocity graph can be converted to a transition matrix by applying a Gaussian kernel to transform the cosine correlations into actual transition probabilities. You can access the Markov transition matrix via `scv.utils.get_transition_matrix`. 
#
# As mentioned, it is internally used to project the velocities into a low-dimensional embedding by applying the mean transition with respect to the transition probabilities, obtained with `scv.tl.velocity_embedding`. Further, we can trace cells along the Markov chain to their origins and potential fates, thereby getting root cells and end points within a trajectory, obtained via `scv.tl.terminal_states`.

# %% [markdown] {"hideCode": false, "hidePrompt": false}
# ### Project the velocities

# %% [markdown] {"hideCode": false, "hidePrompt": false}
# Finally, the velocities are projected onto any embedding, specified by `basis`, and visualized in one of these ways: 
# - on cellular level with `scv.pl.velocity_embedding`,
# - as gridlines with `scv.pl.velocity_embedding_grid`,
# - or as streamlines with `scv.pl.velocity_embedding_stream`.
#
# Note, that the data has an already pre-computed UMAP embedding, and annotated clusters. When applying to your own data, these can be obtained with `scv.tl.umap` and `scv.tl.louvain`. For more details, see the [scanpy tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html). Further, all plotting functions are defaulted to using `basis='umap'` and `color='clusters'`, which you can set accordingly.

# %% {"hideCode": false, "hidePrompt": false}
scv.pl.velocity_embedding_stream(adata, basis='umap')


# %% [markdown] {"hideCode": false, "hidePrompt": false}
# The velocity vector field displayed as streamlines yields fine-grained insights into the developmental processes. It accurately delineates the cycling population of ductal cells and endocrine progenitors. Further, it illuminates cell states of lineage commitment, cell-cycle exit, and endocrine cell differentiation. 

# %% [markdown] {"hideCode": false, "hidePrompt": false}
# The most fine-grained resolution of the velocity vector field we get at single-cell level, with each arrow showing the direction and speed of movement of an individual cell. That reveals, e.g., the early endocrine commitment of Ngn3-cells (yellow) and a clear-cut difference between near-terminal α-cells (blue) and transient β-cells (green).

# %% {"hideCode": false, "hidePrompt": false}
scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120)


# %% [markdown] {"hideCode": false, "hidePrompt": false}
# ### Interpret the velocities

# %% [markdown] {"hideCode": false, "hidePrompt": false}
# This is perhaps the most important part as we advise the user not to limit biological conclusions to the projected velocities, but to examine individual gene dynamics via phase portraits to understand how inferred directions are supported by particular genes. 
#
# See the gif [here](https://user-images.githubusercontent.com/31883718/80227452-eb822480-864d-11ea-9399-56886c5e2785.gif) to get an idea of how to interpret a spliced vs. unspliced phase portrait. Gene activity is orchestrated by transcriptional regulation. Transcriptional induction for a particular gene results in an increase of (newly transcribed) precursor unspliced mRNAs while, conversely, repression or absence of transcription results in a decrease of unspliced mRNAs. Spliced mRNAs is produced from unspliced mRNA and follows the same trend with a time lag. Time is a hidden/latent variable. Thus, the dynamics needs to be inferred from what is actually measured: spliced and unspliced mRNAs as displayed in the phase portrait. 
#
# Now, let us examine the phase portraits of some marker genes, visualized with <br/> 
# `scv.pl.velocity(adata, gene_names)` or `scv.pl.scatter(adata, gene_names)`. 

# %% {"hideCode": false, "hidePrompt": false}
scv.pl.velocity(adata, ['Cpe',  'Gnao1', 'Ins2', 'Adk'], ncols=2)


# %% [markdown] {"hideCode": false, "hidePrompt": false}
# The black line corresponds to the estimated 'steady-state' ratio, i.e. the ratio of unspliced to spliced mRNA abundance which is in a constant transcriptional state. RNA velocity for a particular gene is determined as the residual, i.e. how much an observation deviates from that steady-state line. Positive velocity indicates that a gene is up-regulated, which occurs for cells that show higher abundance of unspliced mRNA for that gene than expected in steady state. Conversely, negative velocity indicates that a gene is down-regulated.
#
# For instance *Cpe* explains the directionality in the up-regulated Ngn3 (yellow) to Pre-endocrine (orange) to β-cells (green), while *Adk*  explains the directionality in the down-regulated Ductal (dark green) to Ngn3 (yellow) to the remaining endocrine cells.

# %% {"hideCode": false, "hidePrompt": false}
scv.pl.scatter(adata, 'Cpe', color=['clusters', 'velocity'], 
               add_outline='Ngn3 high EP, Pre-endocrine, Beta')

# %% [markdown]
# ## Postprocessing

# %% [markdown] {"hideCode": false, "hidePrompt": false}
# ### Identify important genes
# We need a systematic way to identify genes that may help explain the resulting vector field and inferred lineages. 
# To do so, we can test which genes have cluster-specific differential velocity expression, being siginificantly higher/lower compared to the remaining population. The module `scv.tl.rank_velocity_genes` runs a differential velocity t-test and outpus a gene ranking for each cluster. Thresholds can be set (e.g. `min_corr`) to restrict the test on a selection of gene candidates.

# %% {"hideCode": false, "hidePrompt": false}
scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

# %% {"hideCode": false, "hidePrompt": false}
kwargs = dict(frameon=False, size=10, linewidth=1.5, 
              add_outline='Ngn3 high EP, Pre-endocrine, Beta')

scv.pl.scatter(adata, df['Ngn3 high EP'][:5], ylabel='Ngn3 high EP', **kwargs)
scv.pl.scatter(adata, df['Pre-endocrine'][:5], ylabel='Pre-endocrine', **kwargs)


# %% [markdown] {"hideCode": false, "hidePrompt": false}
# The genes *Ptprs, Pclo, Pam, Abcc8, Gnas*, for instance, support the directionality from **Ngn3 high EP** (yellow) to **Pre-endocrine** (orange) to **Beta** (green).

# %% [markdown] {"hideCode": false, "hidePrompt": false}
# ### Velocities in cycling progenitors

# %% [markdown] {"hideCode": false, "hidePrompt": false}
# The cell cycle detected by RNA velocity, is biologically affirmed by cell cycle scores (standardized scores of mean expression levels of phase marker genes).

# %% {"hideCode": false, "hidePrompt": false}
scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])


# %% [markdown] {"hideCode": false, "hidePrompt": false}
# For the cycling Ductal cells, we may screen through S and G2M phase markers. The previous module also computed a spearmans correlation score, which we can use to rank/sort the phase marker genes to then display their phase portraits. 

# %% {"hideCode": false, "hidePrompt": false}
s_genes, g2m_genes = scv.utils.get_phase_marker_genes(adata)
s_genes = scv.get_df(adata[:, s_genes], 'spearmans_score', sort_values=True).index
g2m_genes = scv.get_df(adata[:, g2m_genes], 'spearmans_score', sort_values=True).index

kwargs = dict(frameon=False, ylabel='cell cycle genes')
scv.pl.scatter(adata, list(s_genes[:2]) + list(g2m_genes[:3]), **kwargs)


# %% [markdown] {"hideCode": false, "hidePrompt": false}
# Particularly *Hells* and *Top2a* are well-suited to explain the vector field in the cycling progenitors. 
# *Top2a* gets assigned a high velocity shortly before it actually peaks in the G2M phase. There, the negative velocity then perfectly matches the immediately following down-regulation. 

# %% {"hideCode": false, "hidePrompt": false}
scv.pl.velocity(adata, ['Hells', 'Top2a'], ncols=2, add_outline=True)

# %% [markdown]
# ### Speed and coherence
# Two more useful stats:<br/> 
# - The speed or rate of differentiation is given by the length of the velocity vector. <br/> 
# - The coherence of the vector field (i.e., how a velocity vector correlates with its neighboring velocities) provides a measure of confidence.

# %%
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])


# %% [markdown]
# These provide insights where cells differentiate at a slower/faster pace, and where the direction is un-/determined. 
#
# On cluster-level, we find that differentiation substantially speeds up after cell cycle exit (Ngn3 low EP), keeping the pace during Beta cell production while slowing down during Alpha cell production. 

# %%
df = adata.obs.groupby('clusters')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)


# %% [markdown] {"hideCode": false, "hidePrompt": false}
# ### Velocity graph and pseudotime
# We can visualize the velocity graph to portray all velocity-inferred cell-to-cell connections/transitions. It can be confined to high-probability transitions by setting a `threshold`. The graph, for instance, indicates two phases of Epsilon cell production, coming from early and late Pre-endocrine cells.

# %% {"hideCode": false, "hidePrompt": false}
scv.pl.velocity_graph(adata, threshold=.1)

# %% [markdown]
# Further, the graph can be used to draw descendents/anscestors coming from a specified cell. Here, a pre-endocrine cell is traced to its potential fate.

# %% {"hideCode": false, "hidePrompt": false}
x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)


# %% [markdown]
# Finally, based on the velocity graph, a *velocity pseudotime* can be computed. After inferring a distribution over root cells from the graph, it measures the average number of steps it takes to reach a cell after walking along the graph starting from the root cells. 
#
# Contrarily to diffusion pseudotime, it implicitly infers the root cells and is based on the directed velocity graph instead of the similarity-based diffusion kernel.

# %% {"hideCode": false, "hidePrompt": false}
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')


# %% [markdown]
# ### PAGA velocity graph
# [PAGA](https://doi.org/10.1186/s13059-019-1663-x) graph abstraction has benchmarked as top-performing method for trajectory inference. It provides a graph-like map of the data topology with weighted edges corresponding to the connectivity between two clusters. Here, PAGA is extended by velocity-inferred directionality.

# %%
# PAGA requires to install igraph, if not done yet.
# # !pip install python-igraph --upgrade --quiet

# %% {"hideCode": false, "hidePrompt": false, "tags": []}
# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='clusters')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# %% [markdown]
# This reads from left/row to right/column, thus e.g. assigning a confident transition from Ductal to Ngn3 low EP. 
#
# This table can be summarized by a directed graph superimposed onto the UMAP embedding.

# %%
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)
