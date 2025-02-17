# Advanced Single-Cell RNA-Seq Analysis Pipeline

![Nextflow](https://img.shields.io/badge/nextflow-%23E34F26.svg?style=flat&logo=nextflow&logoColor=white)
![License](https://img.shields.io/badge/license-MIT-blue)

**Unlock the full potential of single-cell RNA sequencing data with a robust, reproducible, and scalable analysis workflow.**  

Single-cell RNA sequencing (scRNA-seq) has revolutionized our ability to study cellular heterogeneity, gene regulation, and developmental trajectories at unprecedented resolution. However, analyzing these complex datasets requires rigorous quality control, sophisticated statistical methods, and computational reproducibility. This pipeline addresses these challenges head-on by integrating **state-of-the-art tools** and **best practices** for end-to-end scRNA-seq analysis—from raw sequencing reads to annotated cell clusters.

## Why This Pipeline?
- **Rigorously Validated QC**: Automated detection of low-quality cells, ambient RNA, and doublets using tools like `SoupX` and `Scrublet`.
- **Reproducible & Scalable**: Built with Nextflow and Docker/Singularity for seamless deployment across laptops, HPCs, and cloud platforms.
- **Advanced Analytics**: Leverages cutting-edge methods like `SCTransform` for normalization, `Harmony` for batch correction, and `SingleR` for cell-type annotation.
- **Benchmark-Ready**: Compare results against public references (e.g., Human Cell Atlas) and validate clusters biologically.

**Key Workflow Steps**:  
1. Raw read QC and trimming  
2. Alignment/quantification (STAR/Alevin)  
3. Cell/gene filtering (adaptive MAD thresholds)  
4. Doublet detection and ambient RNA correction  
5. Normalization, clustering, and annotation  
6. Interactive HTML reports and publication-ready plots  

**Designed For**:  
- Researchers analyzing 10x Genomics, Drop-seq, or Smart-seq2 data  
- Bioinformaticians seeking reproducibility in collaborative projects  
- Core facilities needing standardized workflows  

**Trusted By**:  
[Optional: Add institutional logos or quotes from beta testers]  

---

### 2. Build Singularity Containers (Advanced Tools)
```bash
# Build SCENIC container
singularity build scenic.sif docker://aertslab/pyscenic:0.12.1
```
## Pipeline Steps (Advanced Workflow)

### 1. Ultra-Fast Preprocessing
**Tool**: `kallisto | bustools` (kb-python)  
**Code**:
```bash
kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa GRCh38.fa GRCh38.gtf
kb count -i index.idx -g t2g.txt -x 10xv3 -o out --lamanno --h5ad
```

### 2. Multi-Sample Integration
**Tool**: `Harmony`, `BBKNN`  
**Code** (Seurat v5):
```r
obj <- IntegrateLayers(obj, method = HarmonyIntegration, 
                      orig.reduction = "pca", new.reduction = "harmony")
```

### 3. RNA Velocity Analysis
**Tool**: `scVelo` (Dynamical Model)  
**Code**:
```python
import scvelo as scv
adata = scv.read("spliced_unspliced.h5ad")
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
```

### 4. Regulatory Network Inference (SCENIC)
**Code**:
```bash
# Run in Singularity
singularity exec scenic.sif pyscenic grn \
  --num_workers 20 \
  -o adj.tsv \
  filtered_counts.loom \
  allTFs_hg38.txt
```

### 5. Cell-Cell Communication (NicheNet)
**Code**:
```r
library(nichenetr)
ligand_target_matrix <- readRDS("ligand_target_matrix.rds")
weighted_networks <- readRDS("weighted_networks.rds")
nichenet_output <- predict_ligand_activity(seurat_obj, ligand_target_matrix)
```

## Advanced Analysis Modules

### 1. Multi-Omic Integration (CITE-seq)
```r
# Seurat v5 with TotalVI
library(Seurat)
obj <- CreateSeuratObject(counts = rna_counts, assay = "RNA")
obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)
obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR")
obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "apca"))
```

### 2. Spatial Mapping (Seurat RPCA)
```r
# Reference: scRNA-seq, Query: Visium
anchors <- FindTransferAnchors(
  reference = reference_obj,
  query = spatial_obj,
  normalization.method = "SCT",
  reduction = "rpca"
)
spatial_obj <- MapQuery(
  anchorset = anchors,
  query = spatial_obj,
  reference = reference_obj
)
```

### 3. Machine Learning Annotation (scANVI)
```python
# scANVI semi-supervised training
from scvi.model import SCANVI
scanvi_model = SCANVI.from_rna_model(
  rna_model,
  unlabeled_category="Unknown",
  labels_key="celltype"
)
scanvi_model.train(max_epochs=100)
```

## Usage (Advanced Mode)

### Run with RNA Velocity
```bash
nextflow run main.nf \
  --reads "data/*_R{1,2}.fastq.gz" \
  --modality velocity \
  --splice_info splice_annot.gtf \
  --outdir results_advanced
```

### Run SCENIC + CellChat
```bash
nextflow run main.nf \
  -profile scenic,cellchat \
  --tf_db cisTarget_hg38.feather \
  --interaction_db cellchatdb.rds
```

## Output Structure (Advanced)
```
results_advanced/
├── velocity/            # RNA velocity plots (stream plots)
├── scenic/              # AUCell matrices, TF modules
├── cellchat/            # Communication probability networks
├── spatial_mapping/     # Seurat RPCA transfer results
└── multiomic/           # CITE-seq WNN umap plots
```

## Configuration

### 1. GPU Acceleration (`conf/gpu.config`)
```nextflow
process {
  withLabel 'gpu' {
    clusterOptions = '--gpus=1'
    container = 'nvcr.io/nvidia/pytorch:22.02-py3'
  }
}
```

### 2. Custom Databases
```bash
# Download pre-built SCENIC databases
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
```

## References
1. [scVelo Dynamical Modeling](https://scvelo.readthedocs.io/en/stable/DynamicalModeling/)
2. [SCENIC Protocol](https://www.nature.com/articles/s41596-020-0336-2)
3. [NicheNet Cell Communication](https://doi.org/10.1038/s41592-019-0667-5)

---

### Key Advancements in This Version:
1. **Multi-Modal Integration**: Seurat v5 Weighted Nearest Neighbors (WNN) for CITE-seq
2. **Dynamical RNA Velocity**: Goes beyond steady-state assumptions
3. **Production-Grade SCENIC**: Uses Singularity containers for reproducibility
4. **Semi-Supervised Learning**: scANVI for cell annotation with partial labels
5. **Spatial Mapping**: Integrates scRNA-seq with Visium/ST data
6. **GPU Optimization**: Configurable for NVIDIA GPUs

## Author

**Etienne Ntumba Kabongo**  
📧 Email: [etienne.ntumba.kabongo@umontreal.ca](mailto:etienne.ntumba.kabongo@umontreal.ca)  
🔗 GitHub: [EtienneNtumba](https://github.com/EtienneNtumba)



