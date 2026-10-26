#!/usr/bin/env python3
"""
Step 02: Alpha Diversity and Rarefaction Analysis
Generates:
- Fig2_Alpha_Diversity_and_Rarefaction (.png, .pdf)
- TableS2_Alpha_Diversity_Metrics (.tsv, .csv)
"""

import os
import sys
import zipfile
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from common import (
    set_publication_style, save_figure, save_table,
    load_metadata, SAMPLE_PALETTE, SAMPLE_MARKERS, SAMPLE_ORDER
)

def parse_rarefaction_csv(zip_ref, filename):
    """Parse iterations across depths from alpha-rarefaction.qzv."""
    target = None
    for f in zip_ref.namelist():
        if f.endswith(filename):
            target = f
            break
    if not target:
        return None
    
    df = pd.read_csv(zip_ref.open(target))
    df = df.set_index('sample-id')
    
    # Column format: depth-X_iter-Y
    cols = [c for c in df.columns if 'depth-' in c]
    depths = sorted(list(set(int(re.search(r'depth-(\d+)_', c).group(1)) for c in cols)))
    
    # Compute mean and std for each sample at each depth
    records = []
    for sample in df.index:
        for d in depths:
            d_cols = [c for c in cols if f'depth-{d}_' in c]
            vals = df.loc[sample, d_cols].dropna().astype(float)
            if len(vals) > 0:
                records.append({
                    'sample_id': sample,
                    'depth': d,
                    'mean': vals.mean(),
                    'std': vals.std(),
                    'n': len(vals)
                })
    return pd.DataFrame(records)

def run():
    print("\n--- Running Step 02: Alpha Diversity Analysis ---")
    set_publication_style()

    metadata = load_metadata("metadata.tsv")
    
    import qiime2
    # Prioritize organelle-filtered core metrics (depth 64 retains all 4 samples)
    core_dir = "core-metrics-results-depth64" if os.path.exists("core-metrics-results-depth64") else "core-metrics-results"
    rarefaction_depth = 64 if "depth64" in core_dir else 1413
    print(f"  Using alpha diversity metrics from: {core_dir} (Rarefaction Depth = {rarefaction_depth})")

    metrics_map = {
        'shannon_entropy': f'{core_dir}/shannon_vector.qza',
        'faith_pd': f'{core_dir}/faith_pd_vector.qza',
        'observed_features': f'{core_dir}/observed_features_vector.qza',
        'pielou_evenness': f'{core_dir}/evenness_vector.qza'
    }
    
    dfs = []
    for col_name, path in metrics_map.items():
        art = qiime2.Artifact.load(path)
        df_vec = art.view(qiime2.Metadata).to_dataframe()
        # Find numeric column
        num_col = [c for c in df_vec.columns if c in [col_name, col_name.lower(), 'shannon_entropy', 'faith_pd', 'observed_features', 'pielou_evenness']][0]
        df_vec = df_vec[[num_col]].rename(columns={num_col: col_name})
        dfs.append(df_vec)
        
    alpha_df = pd.concat(dfs, axis=1)
    alpha_df.index.name = 'sample-id'
    
    # Merge with metadata
    full_alpha = metadata.join(alpha_df, how='inner')
    full_alpha = full_alpha.reindex([s for s in SAMPLE_ORDER if s in full_alpha.index])
    
    # Save Table S2
    save_table(full_alpha, "TableS2_Alpha_Diversity_Metrics")
    
    # Load rarefaction curves (prioritize organelle-filtered)
    rare_shannon = None
    rare_obs = None
    rare_qzv = "alpha-rarefaction-no-organelles.qzv" if os.path.exists("alpha-rarefaction-no-organelles.qzv") else "alpha-rarefaction.qzv"
    if os.path.exists(rare_qzv):
        print(f"  Using rarefaction curves from: {rare_qzv}")
        with zipfile.ZipFile(rare_qzv) as z:
            rare_shannon = parse_rarefaction_csv(z, "shannon.csv")
            rare_obs = parse_rarefaction_csv(z, "observed_features.csv")
            
    # Figure 2: Alpha Diversity Box/Bar Panels + Rarefaction Saturation Curves
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.1, 1.2], hspace=0.35, wspace=0.32)
    
    plot_specs = [
        ('shannon_entropy', 'A. Shannon Entropy\n(Richness & Evenness)', (0, 0), 'Entropy (nats)'),
        ('observed_features', 'B. Observed Features\n(ASV Richness)', (0, 1), 'Number of ASVs'),
        ('faith_pd', "C. Faith's Phylogenetic Diversity\n(Branch Length)", (0, 2), "Faith's PD"),
        ('pielou_evenness', "D. Pielou's Evenness\n(Community Uniformity)", (0, 3), "Pielou's J")
    ]
    
    for metric_col, title, grid_pos, y_label in plot_specs:
        ax = fig.add_subplot(gs[grid_pos[0], grid_pos[1]])
        x_pos = np.arange(len(full_alpha))
        vals = full_alpha[metric_col]
        
        bars = ax.bar(x_pos, vals, color=[SAMPLE_PALETTE.get(s, '#333333') for s in full_alpha.index],
                      edgecolor='black', linewidth=0.7, width=0.55, alpha=0.9)
        
        # Add exact value labels on top of bars
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + (max(vals)*0.02),
                    f"{val:.2f}" if metric_col != 'observed_features' else f"{int(val)}",
                    ha='center', va='bottom', fontsize=8, fontweight='bold')
            
        ax.set_title(title, loc='center', fontweight='bold', fontsize=9.5)
        ax.set_ylabel(y_label, fontweight='bold')
        ax.set_xticks(x_pos)
        ax.set_xticklabels([s.replace('_', '\n') for s in full_alpha.index], fontsize=8)
        ax.set_ylim(0, max(vals) * 1.18)
        ax.grid(axis='y', linestyle='--', alpha=0.3)
        
    # Bottom Panels: Rarefaction Curves
    if rare_shannon is not None and rare_obs is not None:
        # Panel E: Shannon Rarefaction
        ax_rare_sh = fig.add_subplot(gs[1, :2])
        for sample in SAMPLE_ORDER:
            s_data = rare_shannon[rare_shannon['sample_id'] == sample].sort_values('depth')
            if not s_data.empty:
                color = SAMPLE_PALETTE.get(sample, '#333333')
                marker = SAMPLE_MARKERS.get(sample, 'o')
                display_name = sample.replace('_', ' ')
                ax_rare_sh.plot(s_data['depth'], s_data['mean'], label=display_name,
                                color=color, marker=marker, markersize=4, linewidth=1.5, alpha=0.9)
                ax_rare_sh.fill_between(s_data['depth'], s_data['mean'] - s_data['std'],
                                        s_data['mean'] + s_data['std'], color=color, alpha=0.15)
                
        ax_rare_sh.set_title("E. Shannon Diversity Rarefaction Saturation Curves", loc='left', fontweight='bold')
        ax_rare_sh.set_xlabel("Sequencing Depth (Reads per Sample)", fontweight='bold')
        ax_rare_sh.set_ylabel("Shannon Entropy", fontweight='bold')
        ax_rare_sh.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
        ax_rare_sh.axvline(rarefaction_depth, color='grey', linestyle=':', linewidth=1.2, label=f'Rarefaction Depth ({rarefaction_depth:,})')
        ax_rare_sh.legend(frameon=True, facecolor='white', framealpha=0.9, edgecolor='none', loc='lower right', fontsize=8)
        ax_rare_sh.grid(True, linestyle='--', alpha=0.3)
        
        # Panel F: Observed Features Rarefaction
        ax_rare_obs = fig.add_subplot(gs[1, 2:])
        for sample in SAMPLE_ORDER:
            s_data = rare_obs[rare_obs['sample_id'] == sample].sort_values('depth')
            if not s_data.empty:
                color = SAMPLE_PALETTE.get(sample, '#333333')
                marker = SAMPLE_MARKERS.get(sample, 'o')
                display_name = sample.replace('_', ' ')
                ax_rare_obs.plot(s_data['depth'], s_data['mean'], label=display_name,
                                 color=color, marker=marker, markersize=4, linewidth=1.5, alpha=0.9)
                ax_rare_obs.fill_between(s_data['depth'], s_data['mean'] - s_data['std'],
                                         s_data['mean'] + s_data['std'], color=color, alpha=0.15)
                
        ax_rare_obs.set_title("F. Observed ASV Richness Rarefaction Curves", loc='left', fontweight='bold')
        ax_rare_obs.set_xlabel("Sequencing Depth (Reads per Sample)", fontweight='bold')
        ax_rare_obs.set_ylabel("Number of Observed ASVs", fontweight='bold')
        ax_rare_obs.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
        ax_rare_obs.axvline(rarefaction_depth, color='grey', linestyle=':', linewidth=1.2, label=f'Rarefaction Depth ({rarefaction_depth:,})')
        ax_rare_obs.legend(frameon=True, facecolor='white', framealpha=0.9, edgecolor='none', loc='lower right', fontsize=8)
        ax_rare_obs.grid(True, linestyle='--', alpha=0.3)
        
    save_figure(fig, "Fig2_Alpha_Diversity_and_Rarefaction")
    plt.close(fig)
    print("Step 02 completed successfully.")

if __name__ == "__main__":
    run()
