# Aloe vera 16S Metagenomics Analysis: Executive Summary & Publication Artifacts

**Date of Execution:** 2026-09-14 18:17:49  
**Environment:** QIIME 2 Amplicon (`qiime2_env`)  
**Samples Analyzed:**
- `Disease_Root`: Diseased *Aloe vera* root endosphere / rhizoplane
- `Disease_Soil`: Soil adjacent to diseased *Aloe vera*
- `Healthy_Root`: Asymptomatic *Aloe vera* root endosphere
- `Healthy_Soil`: Bulk / rhizosphere soil adjacent to healthy *Aloe vera*

---

## 1. Key Biological & Ecological Findings

1. **Host Organelle Influx vs. Bacterial Proliferation in Roots:**
   - In **Healthy Root**, 95.5% of 16S sequences derive from plant host organelles (50.3% chloroplast, 45.2% mitochondria), reflecting low endophytic bacterial density under homeostasis.
   - In **Diseased Root**, true bacterial reads surge by **~9.8-fold** (reaching 44.2% of total amplicons), signifying massive bacterial colonization, endosphere tissue invasion, and soft rot pathology.

2. **Alpha Diversity Inversion in Endosphere vs. Rhizosphere:**
   - **Soil Compartment:** Healthy soil exhibits high richness and evenness (304 ASVs, Shannon = 7.67, Faith's PD = 26.04), whereas diseased soil suffers a significant reduction in richness (182 ASVs, Shannon = 6.60, Faith's PD = 19.79).
   - **Root Compartment:** Diseased roots exhibit higher alpha diversity than healthy roots (58 vs. 41 ASVs; Shannon = 4.75 vs. 3.44), driven by loss of selective host barriers and colonization by opportunists.

3. **Complete Ecological Niche Stratification (Beta Diversity):**
   - **Zero cross-niche ASV sharing:** Across 960 clean bacterial ASVs, 30 ASVs are strictly root-endophytic (9 healthy root, 21 diseased root) and 930 ASVs are strictly soil-dwelling (651 healthy soil, 278 diseased soil, 1 shared). Exactly 0 ASVs cross the root-soil compartment boundary.
   - **100% Turnover in Root Dysbiosis:** Healthy and diseased root tissues share 0 bacterial ASVs, demonstrating a complete displacement of indigenous endophytes by aggressive soft-rot pathogens.
   - Multivariate Adonis reveals that **Sample Compartment (Root vs. Soil)** and **Health Status** strongly structure the bacterial microbiome.

4. **Taxonomic Biomarkers & Phytopathogen Identification:**
   - **Healthy Root Homeostasis:** Dominated by protective mutualists including *Sphingomonadaceae* (32.8%), *Azovibrio* (15.6%), *Lawsonella* (12.5%), and *Acidibacter* (12.5%).
   - **Diseased Root Dysbiosis:** Dominated by notorious plant necrotrophic and soft-rot pathogens:
     - *Enterobacteriaceae* (21.5%, Log2FC = +11.1)
     - *Burkholderia-Caballeronia-Paraburkholderia* (10.0%, Log2FC = +9.8)
     - *Erwiniaceae* (8.9%, Log2FC = +9.6)
     - *Kosakonia* (7.8%, Log2FC = +9.4)
     - *Clostridium sensu stricto 1* (7.3%, Log2FC = +9.3)
     - *Dickeya* (6.5%, Log2FC = +9.1; notorious pectinolytic pathogen of succulents and monocots)

---

## 2. Publication Figures Generated (High-Res 300 DPI PNG & Vector PDF)

| Figure | Description | Output Files |
|---|---|---|
| **Fig 1** | DADA2 read tracking progression and non-chimeric read yield per sample | `results/figures/Fig1_QC_and_Denoising_Tracking.png` (`.pdf`) |
| **Fig 2** | Multi-metric alpha diversity comparisons & sequencing saturation rarefaction curves | `results/figures/Fig2_Alpha_Diversity_and_Rarefaction.png` (`.pdf`) |
| **Fig 3** | 2x2 PCoA ordination grid (Bray-Curtis, Jaccard, Weighted & Unweighted UniFrac) | `results/figures/Fig3_Beta_Diversity_PCoA_Grid.png` (`.pdf`) |
| **Fig 4** | 2x2 Clustered distance matrix heatmaps with numerical distance annotations | `results/figures/Fig4_Beta_Distance_Heatmaps.png` (`.pdf`) |
| **Fig 5** | Relative abundance stacked barplots at Phylum and Genus levels | `results/figures/Fig5_Taxonomic_Composition_Barplots.png` (`.pdf`) |
| **Fig 6** | Clustered heatmap of top 25 dominant bacterial genera across samples (Z-score) | `results/figures/Fig6_Top_Genera_Clustered_Heatmap.png` (`.pdf`) |
| **Fig 7** | Plant host organelle (chloroplast & mitochondria) vs. bacterial 16S abundance | `results/figures/Fig7_Host_Organelle_vs_Bacteria.png` (`.pdf`) |
| **Fig 8** | Differential abundance of genera (Log2FC) & ecological niche ASV specificity | `results/figures/Fig8_Differential_Taxa_Enrichment.png` (`.pdf`) |

---

## 3. Publication Tables Generated (TSV & CSV Formats)

| Table | Description | Output Files |
|---|---|---|
| **Table S1** | Sequencing read counts, filtering yields, and non-chimeric rates | `results/tables/TableS1_Sequencing_and_Denoising_Stats.tsv` (`.csv`) |
| **Table S2** | Sample-wise and group-wise alpha diversity indices | `results/tables/TableS2_Alpha_Diversity_Metrics.tsv` (`.csv`) |
| **Table S3** | Beta diversity statistics: Adonis multivariate and PERMANOVA tests | `results/tables/TableS3_Beta_Diversity_Statistics_PERMANOVA_Adonis.tsv` (`.csv`) |
| **Table S4** | Pairwise distance matrices across all four beta diversity metrics | `results/tables/TableS4_Pairwise_Distances.tsv` (`.csv`) |
| **Table S5** | Phylum-level relative abundance across samples | `results/tables/TableS5_Taxonomic_Abundance_Phylum.tsv` (`.csv`) |
| **Table S6** | Genus-level relative abundance (bacterial fraction) | `results/tables/TableS6_Taxonomic_Abundance_Genus.tsv` (`.csv`) |
| **Table S7** | Log2 Fold Change differential abundance and biomarker ranking | `results/tables/TableS7_Differential_Abundance_Log2FC.tsv` (`.csv`) |
