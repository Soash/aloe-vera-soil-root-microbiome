# Aloe vera 16S Metagenomics Visualization & Analysis Scripts

A modular, publication-ready Python suite for generating publication figures (300 DPI PNG, vector PDF) and structured statistical tables from 16S rRNA amplicon sequencing data.

---

## 1. Environment Setup

The scripts use the QIIME 2 conda environment available on the system:

```bash
# Activate QIIME 2 environment
conda activate qiime2_env
```

Packages utilized:
- `qiime2`
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scikit-bio` (`skbio`)
- `scipy`
- `statsmodels`

---

## 2. Quick Start

Run the entire pipeline with a single command:

```bash
python scripts/run_all.py
```

Or run any step individually:

```bash
python scripts/01_qc_and_denoising.py
python scripts/02_alpha_diversity.py
python scripts/03_beta_diversity.py
python scripts/04_taxonomic_composition.py
python scripts/05_differential_and_biomarkers.py
```

---

## 3. Script Inventory

| Script | Purpose | Key Outputs |
|---|---|---|
| `common.py` | Publication styling (Nature/Cell format, 300 DPI), palettes, artifact loaders | Shared utility module |
| `01_qc_and_denoising.py` | Reads FastQC & DADA2 denoising stats; tracks read progression | `Fig1_QC_and_Denoising_Tracking`, `TableS1_Sequencing_and_Denoising_Stats` |
| `02_alpha_diversity.py` | Alpha diversity indices (Shannon, Faith's PD, Richness, Evenness) & Rarefaction saturation curves | `Fig2_Alpha_Diversity_and_Rarefaction`, `TableS2_Alpha_Diversity_Metrics` |
| `03_beta_diversity.py` | PCoA ordinations (Bray-Curtis, Jaccard, Weighted & Unweighted UniFrac), distance heatmaps, and PERMANOVA/Adonis tests | `Fig3_Beta_Diversity_PCoA_Grid`, `Fig4_Beta_Distance_Heatmaps`, `TableS3_Beta_Diversity_Statistics_PERMANOVA_Adonis`, `TableS4_Pairwise_Distances` |
| `04_taxonomic_composition.py` | Phylum & Genus relative abundance stacked barplots, organelle (chloroplast/mitochondria) vs. bacterial fraction, top 25 genera clustered heatmap | `Fig5_Taxonomic_Composition_Barplots`, `Fig6_Top_Genera_Clustered_Heatmap`, `Fig7_Host_Organelle_vs_Bacteria`, `TableS5_Taxonomic_Abundance_Phylum`, `TableS6_Taxonomic_Abundance_Genus` |
| `05_differential_and_biomarkers.py` | Differential abundance analysis (Log2FC), biomarker discovery for diseased root/soil, and ASV niche specificity | `Fig8_Differential_Taxa_Enrichment`, `TableS7_Differential_Abundance_Log2FC` |
| `run_all.py` | Master workflow runner executing all analysis scripts and compiling `results/summary_report.md` | Executive report & full batch run |

---

## 4. Output Directory Structure

Generated results are saved in `results/`:

```
results/
├── figures/                              # High-resolution (300 DPI) PNG and vector PDF
│   ├── Fig1_QC_and_Denoising_Tracking.png (.pdf)
│   ├── Fig2_Alpha_Diversity_and_Rarefaction.png (.pdf)
│   ├── Fig3_Beta_Diversity_PCoA_Grid.png (.pdf)
│   ├── Fig4_Beta_Distance_Heatmaps.png (.pdf)
│   ├── Fig5_Taxonomic_Composition_Barplots.png (.pdf)
│   ├── Fig6_Top_Genera_Clustered_Heatmap.png (.pdf)
│   ├── Fig7_Host_Organelle_vs_Bacteria.png (.pdf)
│   └── Fig8_Differential_Taxa_Enrichment.png (.pdf)
├── tables/                               # Publication-ready TSV and CSV tables
│   ├── TableS1_Sequencing_and_Denoising_Stats.tsv (.csv)
│   ├── TableS2_Alpha_Diversity_Metrics.tsv (.csv)
│   ├── TableS3_Beta_Diversity_Statistics_PERMANOVA_Adonis.tsv (.csv)
│   ├── TableS4_Pairwise_Distances.tsv (.csv)
│   ├── TableS5_Taxonomic_Abundance_Phylum.tsv (.csv)
│   ├── TableS6_Taxonomic_Abundance_Genus.tsv (.csv)
│   └── TableS7_Differential_Abundance_Log2FC.tsv (.csv)
└── summary_report.md                     # Executive summary report
```
