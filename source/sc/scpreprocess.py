#!/usr/bin/env python

import warnings

warnings.filterwarnings("ignore")

import os
from datetime import datetime

import anndata as ad
import click
import glob
import numpy as np
import pandas as pd
# import ray
import scanpy as sc
import scvi
from scipy.sparse import csr_matrix

# ray.init(num_cpus=4, num_gpus=1)

# @ray.remote(num_gpus=0.25)
def pp(csv_path, ribo_genes, output_directory):
    """
    preprocess csv file containing cell-by-gene count matrices
    """
    # read data
    adata = sc.read_csv(csv_path).T
    
    #'pathwithoutunderscores/GSM5226574_C51ctr_raw_counts.csv'
    sample_id = csv_path.split("_")[1]
    
    print(f"processing: {sample_id}")
    # filter genes in less than 10 cells
    sc.pp.filter_genes(adata, min_cells=10)
    
    # estimate highly variable genes
    sc.pp.highly_variable_genes(
        adata, n_top_genes=2000, subset=True, flavor="seurat_v3"
    )
    
    # estimate doublets with scvi interface to SOLO
    print(f"starting doublet estimation for {sample_id}")
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    df = solo.predict()
    df["prediction"] = solo.predict(soft=False)
    df.index = df.index.map(lambda x: x[:-2])
    df["dif"] = df.doublet - df.singlet
    doublets = df[(df.prediction == "doublet") & (df.dif > 1)]
    print(f"finished doublet estimation for {sample_id}")

    #'pathwithoutunderscores/GSM5226574_C51ctr_raw_counts.csv'
    adata = sc.read_csv(csv_path).T
    adata.obs["Sample"] = sample_id

    # filter doublets
    adata.obs["doublet"] = adata.obs.index.isin(doublets.index)
    adata = adata[~adata.obs.doublet]

    # remove cells with fewer than 200 genes
    sc.pp.filter_cells(adata, min_genes=200)
    # (NOT USED) remove genes that are found in fewer than 3 cells
    # sc.pp.filter_genes( adata, min_cells=3 )

    # annotate the group of mitochondrial genes as 'mt'
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.isin(ribo_genes[0].values)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    )
    
    # filter cells by number of genes, % mitochondrial and % ribosomal
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, 0.98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]
    adata = adata[adata.obs.pct_counts_mt < 20]
    adata = adata[adata.obs.pct_counts_ribo < 2]
    
    # save a copy of the sample dataframe
    adata.write(f"{output_directory}/adata_{sample_id}.h5ad", compression="gzip")
    print(f"wrote {sample_id} output to {output_directory}/adata_{sample_id}.h5ad")
    # TODO: used for testing
    # adata = test_data()
    
    return adata


def test_data():
    counts = csr_matrix(np.random.poisson(1, size=(4, 9)), dtype=np.float32)
    adata = ad.AnnData(counts)
    return adata


@click.command()
@click.option(
    "--input-directory",
    prompt="Enter the directory containing input files",
    help="Directory containing input files of potentially gzipped\n\
          single-cell RNA-seq cell-by-gene count matrices",
)
def cli(input_directory="data/GSE171524/supplementary/"):
    """
    preprocess single-cell RNA-seq data count matrices across samples
    and combine into one file
    """
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")
    # timestamp = "20220930_215746"

    output_directory = f"output/{timestamp}"

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        print(f"\ncreated output directory: {output_directory}")
    else:
        print(f"output directory exists: {output_directory}")

    print(f"input directory: {input_directory}\n")
    
    ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
    ribo_genes = pd.read_table(ribo_url, skiprows=2, header = None)
    ribo_genes

    # filenames = os.listdir(input_directory)
    filenames = [os.path.basename(x) for x in glob.glob(input_directory+"*.csv*")]

    print(filenames)
    
    out = []
    
    # run in serial
    for file in filenames:
        sample_id = file.split("_")[1]
        outfile = f"{output_directory}/adata_{sample_id}.h5ad"
        if os.path.exists(outfile):
            print(f"{outfile} already exists")
            out.append(sc.read_h5ad(outfile))
        else:
            print(f"processing {file}")
            out.append(pp(input_directory + file, ribo_genes, output_directory))
    
    # # run in parallel with ray
    # def parmap(f, list):
    #    return [f.remote(input_directory + x, ribo_genes, output_directory) for x in list]
    
    # result_ids = parmap(pp, filenames)
    # out = ray.get(result_ids)
    
    print("combining AnnData objects")
    adata = sc.concat(out)
    
    print("writing combined AnnData h5ad")
    adata.write(f"{output_directory}/adata_combined.h5ad", compression="gzip")
    print(f"wrote combined data to {output_directory}/adata_combined.h5ad")

if __name__ == "__main__":
    cli()
