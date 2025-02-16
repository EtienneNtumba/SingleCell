# main.py
# Single-Cell RNA-seq Analysis using Scanpy

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import yaml
import os

# Load configuration
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

DATASET_PATH = config['data']['dataset_path']
FILTER_CELLS = config['qc']['min_cells']
FILTER_GENES = config['qc']['min_genes']
NORMALIZATION_METHOD = config['preprocessing']['normalization']
N_PCS = config['dim_reduction']['n_pcs']
CLUSTER_RES = config['clustering']['resolution']
SAVE_PATH = config['output']['save_dir']

def preprocess_data():
    """Load and preprocess scRNA-seq data."""
    print("Loading data...")
    adata = sc.read_10x_mtx(DATASET_PATH, var_names='gene_symbols', cache=True)
    
    print("Performing Quality Control...")
    sc.pp.filter_cells(adata, min_genes=FILTER_CELLS)
    sc.pp.filter_genes(adata, min_cells=FILTER_GENES)
    
    print("Normalizing data...")
    if NORMALIZATION_METHOD == 'log1p':
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif NORMALIZATION_METHOD == 'sctransform':
        sc.experimental.pp.sctransform(adata)
    
    print("Identifying highly variable genes...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    
    return adata

def perform_dimensionality_reduction(adata):
    """Reduce dimensionality using PCA and UMAP."""
    print("Running PCA...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    
    print("Computing nearest neighbors...")
    sc.pp.neighbors(adata, n_pcs=N_PCS)
    
    print("Performing UMAP...")
    sc.tl.umap(adata)
    return adata

def cluster_cells(adata):
    """Cluster single cells into distinct populations."""
    print("Clustering cells...")
    sc.tl.leiden(adata, resolution=CLUSTER_RES)
    return adata

def run_marker_analysis(adata):
    """Identify differentially expressed genes."""
    print("Running differential expression analysis...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    return adata

def visualize_results(adata):
    """Generate visualizations for clustering and expression."""
    print("Generating plots...")
    sc.pl.umap(adata, color=['leiden'], save='_clusters.png')
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, save='_DEGs.png')

def main():
    """Run the full pipeline."""
    adata = preprocess_data()
    adata = perform_dimensionality_reduction(adata)
    adata = cluster_cells(adata)
    adata = run_marker_analysis(adata)
    visualize_results(adata)
    adata.write(os.path.join(SAVE_PATH, "processed_data.h5ad"))
    print("Pipeline completed successfully.")

if __name__ == "__main__":
    main()
