#!/usr/bin/env python3
"""
Step 03: Beta Diversity Ordination and Distance Matrix Analysis
Generates:
- Fig3_Beta_Diversity_PCoA_Grid (.png, .pdf)
- Fig4_Beta_Distance_Heatmaps (.png, .pdf)
- TableS3_Beta_Diversity_Statistics_PERMANOVA_Adonis (.tsv, .csv)
- TableS4_Pairwise_Distances (.tsv, .csv)
"""

import os
import sys
import io
import zipfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from common import (
    set_publication_style, save_figure, save_table,
    load_metadata, SAMPLE_PALETTE, SAMPLE_MARKERS, SAMPLE_ORDER
)

def extract_permanova_stats(qzv_path):
    """Extract test statistic and p-value from PERMANOVA/PERMDISP QZV index.html."""
    with zipfile.ZipFile(qzv_path) as z:
        for f in z.namelist():
            if f.endswith('index.html'):
                html = z.open(f).read().decode('utf-8')
                dfs = pd.read_html(io.StringIO(html))
                t0 = dfs[0].set_index(dfs[0].columns[0])
                stat = float(t0.loc['test statistic'].iloc[0])
                pval = float(t0.loc['p-value'].iloc[0])
                method = str(t0.loc['method name'].iloc[0])
                stat_name = str(t0.loc['test statistic name'].iloc[0])
                return {'method': method, 'statistic_name': stat_name, 'statistic': stat, 'p_value': pval}
    return None

def extract_adonis_table(qzv_path):
    """Extract Adonis regression table from QZV."""
    with zipfile.ZipFile(qzv_path) as z:
        for f in z.namelist():
            if f.endswith('adonis.tsv'):
                df = pd.read_csv(z.open(f), sep='\t', index_col=0)
                return df
    return None

def run():
    print("\n--- Running Step 03: Beta Diversity Analysis ---")
    set_publication_style()

    metadata = load_metadata("metadata.tsv")
    import qiime2
    from skbio import OrdinationResults, DistanceMatrix

    core_dir = "core-metrics-results-depth64" if os.path.exists("core-metrics-results-depth64") else "core-metrics-results"
    print(f"  Using beta diversity metrics from: {core_dir}")

    metrics = [
        ('bray_curtis', 'Bray-Curtis (Abundance-based Dissimilarity)'),
        ('jaccard', 'Jaccard (Incidence-based Dissimilarity)'),
        ('weighted_unifrac', 'Weighted UniFrac (Quantitative Phylogenetic)'),
        ('unweighted_unifrac', 'Unweighted UniFrac (Qualitative Phylogenetic)')
    ]

    pcoa_data = {}
    dist_matrices = {}

    for m_key, m_label in metrics:
        pcoa_file = f"{core_dir}/{m_key}_pcoa_results.qza"
        dm_file = f"{core_dir}/{m_key}_distance_matrix.qza"

        art_pcoa = qiime2.Artifact.load(pcoa_file)
        ord_res = art_pcoa.view(OrdinationResults)
        pcoa_data[m_key] = ord_res

        art_dm = qiime2.Artifact.load(dm_file)
        dm = art_dm.view(DistanceMatrix)
        dist_matrices[m_key] = dm.to_data_frame()

    # Save Table S4: Pairwise Distances for all metrics
    dm_combined_rows = []
    for m_key, m_label in metrics:
        df_dm = dist_matrices[m_key]
        for s1 in df_dm.index:
            for s2 in df_dm.columns:
                if s1 < s2:  # unique pairs
                    dm_combined_rows.append({
                        'Metric': m_key,
                        'Sample_1': s1,
                        'Sample_2': s2,
                        'Distance': df_dm.loc[s1, s2]
                    })
    save_table(pd.DataFrame(dm_combined_rows), "TableS4_Pairwise_Distances", index=False)

    # Compile Table S3: Statistical Tests (Adonis & PERMANOVA)
    stat_records = []
    for m_key, m_label in metrics:
        # Adonis
        adonis_qzv = f"{core_dir}/{m_key}-adonis.qzv"
        if os.path.exists(adonis_qzv):
            ad_df = extract_adonis_table(adonis_qzv)
            for factor in ['Health_Status', 'Sample_Type', 'Residuals', 'Total']:
                if factor in ad_df.index:
                    stat_records.append({
                        'Metric': m_key,
                        'Test': 'Adonis (Multivariate)',
                        'Factor': factor,
                        'Df': ad_df.loc[factor, 'Df'],
                        'SumsOfSqs': ad_df.loc[factor, 'SumsOfSqs'],
                        'R2': ad_df.loc[factor, 'R2'],
                        'F_statistic': ad_df.loc[factor, 'F.Model'] if 'F.Model' in ad_df.columns else np.nan,
                        'p_value': ad_df.loc[factor, 'Pr(>F)'] if 'Pr(>F)' in ad_df.columns else np.nan
                    })
        # PERMANOVA
        for col in ['Health_Status', 'Sample_Type']:
            perm_qzv = f"{core_dir}/{m_key}-{col}-permanova.qzv"
            if os.path.exists(perm_qzv):
                res = extract_permanova_stats(perm_qzv)
                if res:
                    stat_records.append({
                        'Metric': m_key,
                        'Test': 'PERMANOVA (Univariate)',
                        'Factor': col,
                        'Df': 1,
                        'SumsOfSqs': np.nan,
                        'R2': np.nan,
                        'F_statistic': res['statistic'],
                        'p_value': res['p_value']
                    })
    save_table(pd.DataFrame(stat_records), "TableS3_Beta_Diversity_Statistics_PERMANOVA_Adonis", index=False)

    # Figure 3: 2x2 PCoA Grid
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for i, (m_key, m_label) in enumerate(metrics):
        ax = axes[i]
        ord_res = pcoa_data[m_key]
        coords = ord_res.samples
        prop_exp = ord_res.proportion_explained

        pc1_pct = prop_exp[0] * 100
        pc2_pct = prop_exp[1] * 100

        for sample in SAMPLE_ORDER:
            if sample in coords.index:
                x = coords.loc[sample, 0]
                y = coords.loc[sample, 1]
                color = SAMPLE_PALETTE.get(sample, '#333333')
                marker = SAMPLE_MARKERS.get(sample, 'o')
                ax.scatter(x, y, color=color, marker=marker, s=160, edgecolor='black', linewidth=1.2, zorder=4)
                
                # Annotate point label with anti-collision offset
                label_text = sample.replace('_', ' ')
                offset_xy = (-10, -12) if (sample == 'Disease_Soil' and m_key == 'unweighted_unifrac') else (8, 6)
                ha_val = 'right' if (sample == 'Disease_Soil' and m_key == 'unweighted_unifrac') else 'left'
                ax.annotate(label_text, (x, y), textcoords="offset points",
                            xytext=offset_xy, ha=ha_val, fontsize=8.5, fontweight='bold',
                            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.75))

        ax.axhline(0, color='grey', linestyle='--', linewidth=0.8, alpha=0.5, zorder=1)
        ax.axvline(0, color='grey', linestyle='--', linewidth=0.8, alpha=0.5, zorder=1)
        
        # Expand limits slightly for labels
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        x_margin = max((x_max - x_min) * 0.25, 0.1)
        y_margin = max((y_max - y_min) * 0.25, 0.1)
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)

        ax.set_title(f"{chr(65+i)}. {m_label}", loc='left', fontweight='bold', fontsize=10)
        ax.set_xlabel(f"PCoA 1 ({pc1_pct:.1f}% explained)", fontweight='bold')
        ax.set_ylabel(f"PCoA 2 ({pc2_pct:.1f}% explained)", fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.25)

    plt.tight_layout()
    save_figure(fig, "Fig3_Beta_Diversity_PCoA_Grid")
    plt.close(fig)

    # Figure 4: Beta Distance Heatmaps (2x2)
    fig, axes = plt.subplots(2, 2, figsize=(11, 9.5))
    axes = axes.flatten()

    for i, (m_key, m_label) in enumerate(metrics):
        ax = axes[i]
        df_dm = dist_matrices[m_key].reindex(index=SAMPLE_ORDER, columns=SAMPLE_ORDER)
        labels = [s.replace('_', ' ') for s in SAMPLE_ORDER]

        sns.heatmap(df_dm, ax=ax, annot=True, fmt=".3f", cmap="YlOrRd", cbar_kws={'label': 'Distance'},
                    xticklabels=labels, yticklabels=labels, linewidths=0.5, linecolor='white',
                    annot_kws={"size": 9, "weight": "bold"})

        ax.set_title(f"{chr(65+i)}. {m_label}", loc='left', fontweight='bold', fontsize=9.5)
        ax.set_xticklabels(labels, rotation=30, ha='right', fontweight='bold')
        ax.set_yticklabels(labels, rotation=0, fontweight='bold')

    plt.tight_layout()
    save_figure(fig, "Fig4_Beta_Distance_Heatmaps")
    plt.close(fig)
    print("Step 03 completed successfully.")

if __name__ == "__main__":
    run()
