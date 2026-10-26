#!/usr/bin/env python3
"""
Step 01: Sequencing Quality Control and DADA2 Denoising Yield Analysis
Generates:
- Fig1_QC_and_Denoising_Tracking (.png, .pdf)
- TableS1_Sequencing_and_Denoising_Stats (.tsv, .csv)
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from common import (
    set_publication_style, save_figure, save_table,
    load_metadata, SAMPLE_PALETTE, SAMPLE_ORDER
)

def run():
    print("\n--- Running Step 01: QC and Denoising Analysis ---")
    set_publication_style()

    metadata = load_metadata("metadata.tsv")
    
    # Load denoising stats
    import qiime2
    denois_art = qiime2.Artifact.load("denoising-stats.qza")
    denois_df = denois_art.view(qiime2.Metadata).to_dataframe()
    
    # Cast numeric columns
    for col in denois_df.columns:
        denois_df[col] = pd.to_numeric(denois_df[col], errors='coerce')
    
    # Merge with metadata
    merged_df = metadata.join(denois_df, how='inner')
    merged_df = merged_df.reindex([s for s in SAMPLE_ORDER if s in merged_df.index])
    
    # Check for organelle-filtered table to report clean bacterial yields
    if os.path.exists("table-no-organelles.qza"):
        tbl_clean_art = qiime2.Artifact.load("table-no-organelles.qza")
        tbl_clean = tbl_clean_art.view(pd.DataFrame)
        clean_reads = tbl_clean.sum(axis=1)
        merged_df['bacterial-reads'] = clean_reads
        merged_df['organelle-reads'] = merged_df['non-chimeric'] - clean_reads
        merged_df['organelle-contamination-%'] = ((merged_df['organelle-reads'] / merged_df['non-chimeric']) * 100).round(2)
        merged_df['bacterial-yield-%'] = ((merged_df['bacterial-reads'] / merged_df['input']) * 100).round(2)

    # Save Table S1
    table_s1 = merged_df.copy()
    table_s1.index.name = "Sample_ID"
    save_table(table_s1, "TableS1_Sequencing_and_Denoising_Stats")
    
    # Visualization: Figure 1
    fig = plt.figure(figsize=(12, 5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=0.3)
    
    # Panel A: Read count progression through DADA2 stages
    ax1 = fig.add_subplot(gs[0, 0])
    stages = ['input', 'filtered', 'denoised', 'merged', 'non-chimeric']
    stage_labels = ['Input', 'Filtered', 'Denoised', 'Merged', 'Non-Chimeric']
    
    x = np.arange(len(stage_labels))
    width = 0.18
    offsets = np.linspace(-1.5 * width, 1.5 * width, len(merged_df))
    
    for i, (sample, row) in enumerate(merged_df.iterrows()):
        counts = [row[stage] for stage in stages]
        color = SAMPLE_PALETTE.get(sample, '#333333')
        display_name = sample.replace('_', ' ')
        bars = ax1.bar(x + offsets[i], counts, width, label=display_name, color=color, edgecolor='black', linewidth=0.6, alpha=0.9)
        
    ax1.set_ylabel("Read Pairs Count (Reads)", fontweight='bold')
    ax1.set_title("A. DADA2 Read Processing Pipeline Tracking", loc='left', fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(stage_labels, fontweight='bold')
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{int(y):,}'))
    ax1.legend(frameon=True, facecolor='white', framealpha=0.9, edgecolor='none', loc='upper right')
    ax1.grid(axis='y', linestyle='--', alpha=0.3)
    
    # Panel B: Final non-chimeric reads and retention efficiency
    ax2 = fig.add_subplot(gs[0, 1])
    retention_pct = merged_df['percentage of input non-chimeric']
    non_chimeric_counts = merged_df['non-chimeric']
    
    y_pos = np.arange(len(merged_df))
    bars = ax2.barh(y_pos, retention_pct, color=[SAMPLE_PALETTE.get(s, '#333333') for s in merged_df.index],
                    edgecolor='black', linewidth=0.7, height=0.55, alpha=0.9)
    
    for bar, count, pct in zip(bars, non_chimeric_counts, retention_pct):
        ax2.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                 f"{pct:.2f}% ({int(count):,} reads)",
                 va='center', ha='left', fontsize=8, fontweight='bold', color='#222222')
        
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([s.replace('_', ' ') for s in merged_df.index], fontweight='bold')
    ax2.set_xlabel("Non-Chimeric Read Retention (% of Input)", fontweight='bold')
    ax2.set_xlim(0, max(retention_pct) * 1.5)
    ax2.set_title("B. Overall Non-Chimeric Yield & Retention", loc='left', fontweight='bold')
    ax2.grid(axis='x', linestyle='--', alpha=0.3)
    
    pass
    save_figure(fig, "Fig1_QC_and_Denoising_Tracking")
    plt.close(fig)
    print("Step 01 completed successfully.")

if __name__ == "__main__":
    run()
