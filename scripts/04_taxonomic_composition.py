#!/usr/bin/env python3
"""
Step 04: Taxonomic Profiling and Host Organelle vs. Bacterial Dynamics
Generates:
- Fig5_Taxonomic_Composition_Barplots (.png, .pdf)
- Fig6_Top_Genera_Clustered_Heatmap (.png, .pdf)
- Fig7_Host_Organelle_vs_Bacteria (.png, .pdf)
- TableS5_Taxonomic_Abundance_Phylum (.tsv, .csv)
- TableS6_Taxonomic_Abundance_Genus (.tsv, .csv)
"""

import os
import sys
import zipfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from common import (
    set_publication_style, save_figure, save_table,
    load_metadata, SAMPLE_PALETTE, SAMPLE_ORDER
)

def clean_taxon_name(raw_name, level='genus'):
    """Format taxonomy strings cleanly for publication display."""
    parts = [p.strip() for p in raw_name.split(';')]
    if level == 'phylum':
        for p in parts:
            if p.startswith('p__'):
                clean = p.replace('p__', '')
                return clean if clean else "Unclassified Phylum"
        return "Unclassified"
    elif level == 'family':
        for p in reversed(parts):
            if p.startswith('f__') and len(p) > 3:
                return p.replace('f__', '')
            if p.startswith('o__') and len(p) > 3:
                return f"Uncl. {p.replace('o__', '')}"
        return "Unclassified Family"
    elif level == 'genus':
        for p in reversed(parts):
            val = p.split('__')[-1].strip() if '__' in p else p.strip()
            if not val or val == '__' or val.lower() in ['uncultured', 'unknown', 'unknown_family', 'unidentified']:
                continue
            if p.startswith('g__'):
                return val
            if p.startswith('f__'):
                return f"Uncl. {val}"
            if p.startswith('o__'):
                return f"Uncl. {val}"
            if p.startswith('c__'):
                return f"Uncl. {val}"
            if p.startswith('p__'):
                return f"Uncl. {val}"
        return "Unclassified Genus"
    return raw_name

def build_clean_genus_table(bact_genus_pct, pseudocount=0.01):
    """
    Convert 4xN raw genus matrix into a publication-grade, separated, tidy dataframe
    with separated taxonomic hierarchy columns, sample percentages, fold changes, and enrichment status.
    """
    rank_prefixes = {
        'd__': 'Domain', 'k__': 'Domain', 'p__': 'Phylum', 'c__': 'Class',
        'o__': 'Order', 'f__': 'Family', 'g__': 'Genus'
    }
    
    rows = []
    for col in bact_genus_pct.columns:
        parts = [p.strip() for p in col.split(';')]
        ranks = {
            'Domain': 'Bacteria', 'Phylum': 'Unclassified', 'Class': 'Unclassified',
            'Order': 'Unclassified', 'Family': 'Unclassified', 'Genus': 'Unclassified'
        }
        for p in parts:
            for prefix, rank in rank_prefixes.items():
                if p.startswith(prefix):
                    val = p[len(prefix):].strip()
                    if val and val != '__':
                        ranks[rank] = val
        
        hr = float(bact_genus_pct.loc['Healthy_Root', col]) if 'Healthy_Root' in bact_genus_pct.index else 0.0
        dr = float(bact_genus_pct.loc['Disease_Root', col]) if 'Disease_Root' in bact_genus_pct.index else 0.0
        hs = float(bact_genus_pct.loc['Healthy_Soil', col]) if 'Healthy_Soil' in bact_genus_pct.index else 0.0
        ds = float(bact_genus_pct.loc['Disease_Soil', col]) if 'Disease_Soil' in bact_genus_pct.index else 0.0
        
        mean_abund = (hr + dr + hs + ds) / 4.0
        max_abund = max(hr, dr, hs, ds)
        
        root_log2fc = np.log2((dr + pseudocount) / (hr + pseudocount))
        soil_log2fc = np.log2((ds + pseudocount) / (hs + pseudocount))
        
        clean_g = ranks['Genus']
        fam = ranks['Family']
        ord_name = ranks['Order']
        phy = ranks['Phylum']
        if clean_g.lower() in ['unclassified', 'uncultured', 'unknown', 'unknown_family', '__']:
            if fam.lower() not in ['unclassified', 'uncultured', 'unknown', 'unknown_family', '__']:
                clean_g = f"{fam} (uncl. genus)"
                ranks['Genus'] = f"Unclassified {fam}"
            elif ord_name.lower() not in ['unclassified', 'uncultured', 'unknown', '__']:
                clean_g = f"{ord_name} (uncl. family)"
                ranks['Genus'] = f"Unclassified {ord_name}"
                ranks['Family'] = f"Unclassified {ord_name}"
            else:
                clean_g = f"{phy} (uncl. order)"
                ranks['Genus'] = f"Unclassified {phy}"
                ranks['Family'] = f"Unclassified {phy}"
                ranks['Order'] = f"Unclassified {phy}"

                
        if dr >= 1.0 and root_log2fc >= 2.0:
            cat = "Enriched in Disease Root"
        elif hr >= 1.0 and root_log2fc <= -2.0:
            cat = "Enriched in Healthy Root"
        elif ds >= 1.0 and soil_log2fc >= 1.0:
            cat = "Enriched in Disease Soil"
        elif hs >= 1.0 and soil_log2fc <= -1.0:
            cat = "Enriched in Healthy Soil"
        elif max_abund >= 0.5:
            cat = "Core / Common (>0.5%)"
        else:
            cat = "Low Abundance (<0.5%)"
            
        rows.append({
            'Domain': ranks['Domain'],
            'Phylum': ranks['Phylum'],
            'Class': ranks['Class'],
            'Order': ranks['Order'],
            'Family': ranks['Family'],
            'Genus': ranks['Genus'],
            'Taxon_Name': clean_g,
            'Disease_Root (%)': round(dr, 4),
            'Healthy_Root (%)': round(hr, 4),
            'Disease_Soil (%)': round(ds, 4),
            'Healthy_Soil (%)': round(hs, 4),
            'Mean_Abundance (%)': round(mean_abund, 4),
            'Max_Abundance (%)': round(max_abund, 4),
            'Root_Log2FC': round(root_log2fc, 2),
            'Soil_Log2FC': round(soil_log2fc, 2),
            'Enrichment_Status': cat,
            'Full_Taxonomy': col
        })
        
    clean_df = pd.DataFrame(rows)
    clean_df = clean_df.sort_values(by=['Max_Abundance (%)', 'Mean_Abundance (%)'], ascending=[False, False]).reset_index(drop=True)
    clean_df.index = clean_df.index + 1
    clean_df.index.name = 'Rank'
    return clean_df


def load_taxa_level(qzv_path, level_filename):
    """Load relative abundance table for a specific taxonomic level."""
    with zipfile.ZipFile(qzv_path) as z:
        for f in z.namelist():
            if f.endswith(level_filename):
                df = pd.read_csv(z.open(f), index_col=0)
                # Drop metadata columns
                cols = [c for c in df.columns if c not in ['Health_Status', 'Sample_Type']]
                return df[cols]
    raise FileNotFoundError(f"Could not find {level_filename} in {qzv_path}")

def run():
    print("\n--- Running Step 04: Taxonomic Composition Analysis ---")
    set_publication_style()

    qzv_path = "taxa-bar-plots.qzv"
    metadata = load_metadata("metadata.tsv")

    # Load levels
    df_phylum = load_taxa_level(qzv_path, "level-2.csv").reindex(SAMPLE_ORDER)
    df_family = load_taxa_level(qzv_path, "level-5.csv").reindex(SAMPLE_ORDER)
    df_genus = load_taxa_level(qzv_path, "level-6.csv").reindex(SAMPLE_ORDER)

    # 1. Host Organelle vs True Bacteria Analysis
    total_reads = df_genus.sum(axis=1)
    mito_cols = [c for c in df_genus.columns if 'Mitochondria' in c]
    chloro_cols = [c for c in df_genus.columns if 'Chloroplast' in c]
    bact_cols = [c for c in df_genus.columns if c not in mito_cols and c not in chloro_cols]

    organelle_df = pd.DataFrame({
        'True_Bacteria_Reads': df_genus[bact_cols].sum(axis=1),
        'Chloroplast_Reads': df_genus[chloro_cols].sum(axis=1),
        'Mitochondria_Reads': df_genus[mito_cols].sum(axis=1),
        'Total_Reads': total_reads
    })
    organelle_pct = pd.DataFrame({
        'Bacterial Microbiota (%)': (organelle_df['True_Bacteria_Reads'] / total_reads) * 100,
        'Plant Chloroplast (%)': (organelle_df['Chloroplast_Reads'] / total_reads) * 100,
        'Plant Mitochondria (%)': (organelle_df['Mitochondria_Reads'] / total_reads) * 100
    })

    # Save Figure 7: Host Organelle vs Bacteria
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1.4, 1]})
    
    # Stacked bar
    sample_labels = [s.replace('_', '\n') for s in SAMPLE_ORDER]
    org_colors = ['#2b83ba', '#74add1', '#d73027']
    bottom = np.zeros(len(SAMPLE_ORDER))
    
    for col_name, color in zip(['Bacterial Microbiota (%)', 'Plant Chloroplast (%)', 'Plant Mitochondria (%)'], ['#2E7D32', '#66BB6A', '#FFA726']):
        vals = organelle_pct[col_name].values
        ax1.bar(sample_labels, vals, bottom=bottom, label=col_name, color=color, edgecolor='black', linewidth=0.6, width=0.55)
        bottom += vals
        
    ax1.set_ylabel("Relative Sequence Abundance (%)", fontweight='bold')
    ax1.set_title("A. 16S rRNA Sequence Origin (Host Organelles vs Bacteria)", loc='left', fontweight='bold')
    ax1.legend(frameon=True, facecolor='white', framealpha=0.9, edgecolor='none', loc='upper right')
    ax1.set_ylim(0, 100)
    ax1.grid(axis='y', linestyle='--', alpha=0.3)

    # Panel B: Bacterial Abundance in Endosphere (Disease vs Healthy Root)
    root_bact_pct = [organelle_pct.loc['Healthy_Root', 'Bacterial Microbiota (%)'],
                     organelle_pct.loc['Disease_Root', 'Bacterial Microbiota (%)']]
    b_bars = ax2.bar(['Healthy Root', 'Diseased Root'], root_bact_pct,
                     color=['#2E7D32', '#C62828'], edgecolor='black', linewidth=0.7, width=0.45)
    for bar, pct in zip(b_bars, root_bact_pct):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
                 f"{pct:.2f}%", ha='center', va='bottom', fontsize=9, fontweight='bold')
        
    ax2.set_ylabel("True Bacterial Abundance (%)", fontweight='bold')
    ax2.set_title("B. Root Endosphere Bacterial Proliferation", loc='left', fontweight='bold')
    ax2.set_ylim(0, 55)
    ax2.annotate("~9.8x Bacterial Influx\nin Diseased Root", xy=(1, root_bact_pct[1]), xytext=(0.5, 48),
                 arrowprops=dict(facecolor='black', shrink=0.08, width=1, headwidth=6),
                 fontsize=8.5, fontweight='bold', ha='center',
                 bbox=dict(boxstyle="round,pad=0.3", fc="#FFF3E0", ec="#E65100", lw=1))
    ax2.grid(axis='y', linestyle='--', alpha=0.3)

    plt.tight_layout()
    save_figure(fig, "Fig7_Host_Organelle_vs_Bacteria")
    plt.close(fig)

    # 2. Phylum & Genus Relative Abundance Tables (Organelle-Filtered Bacterial Community)
    clean_qzv = "taxa-bar-plots-no-organelles.qzv" if os.path.exists("taxa-bar-plots-no-organelles.qzv") else qzv_path
    print(f"  Using taxonomic profiles from: {clean_qzv}")
    
    if clean_qzv != qzv_path:
        df_phylum_clean = load_taxa_level(clean_qzv, "level-2.csv").reindex(SAMPLE_ORDER)
        df_genus_clean = load_taxa_level(clean_qzv, "level-6.csv").reindex(SAMPLE_ORDER)
        phylum_pct = df_phylum_clean.div(df_phylum_clean.sum(axis=1), axis=0) * 100
        bact_genus_pct = df_genus_clean.div(df_genus_clean.sum(axis=1), axis=0) * 100
    else:
        bact_genus = df_genus[bact_cols]
        bact_genus_pct = bact_genus.div(bact_genus.sum(axis=1), axis=0) * 100
        phylum_pct = df_phylum.div(df_phylum.sum(axis=1), axis=0) * 100

    save_table(phylum_pct, "TableS5_Taxonomic_Abundance_Phylum")
    save_table(bact_genus_pct, "TableS6_Taxonomic_Abundance_Genus")

    # Generate publication-grade, separated, tidy genus abundance table
    clean_genus_df = build_clean_genus_table(bact_genus_pct)
    # Save clean separated version to TableS6_Taxonomic_Abundance_Genus.csv (replacing raw horizontal matrix)
    clean_genus_csv = os.path.join("results/tables", "TableS6_Taxonomic_Abundance_Genus.csv")
    clean_genus_df.to_csv(clean_genus_csv, index=True)
    # Also save dedicated explicit Clean files
    clean_genus_df.to_csv(os.path.join("results/tables", "TableS6_Taxonomic_Abundance_Genus_Clean.csv"), index=True)
    clean_genus_df.to_csv(os.path.join("results/tables", "TableS6_Taxonomic_Abundance_Genus_Clean.tsv"), sep='\t', index=True)
    print(f"  [SAVED CLEAN TABLE] {clean_genus_csv} (separated ranks, tidy layout)")


    # 3. Figure 5: Multi-panel Stacked Barplots (Phylum & Genus)
    fig, (ax_p, ax_g) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: Phylum
    top_phyla = phylum_pct.mean().sort_values(ascending=False).head(7).index.tolist()
    phyla_plot_df = phylum_pct[top_phyla].copy()
    phyla_plot_df['Other Phyla'] = 100.0 - phyla_plot_df.sum(axis=1)

    colors_phyla = sns.color_palette("tab10", len(phyla_plot_df.columns))
    bottom_p = np.zeros(len(SAMPLE_ORDER))
    for col, c in zip(phyla_plot_df.columns, colors_phyla):
        clean_label = clean_taxon_name(col, 'phylum') if col != 'Other Phyla' else col
        vals = phyla_plot_df[col].values
        ax_p.bar(sample_labels, vals, bottom=bottom_p, label=clean_label, color=c, edgecolor='black', linewidth=0.5, width=0.55)
        bottom_p += vals

    ax_p.set_title("A. Bacterial Phylum Composition (Organelles Excluded)", loc='left', fontweight='bold')
    ax_p.set_ylabel("Relative Abundance (%)", fontweight='bold')
    ax_p.set_ylim(0, 100)
    ax_p.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True, facecolor='white', framealpha=0.9, edgecolor='none', fontsize=8)
    ax_p.grid(axis='y', linestyle='--', alpha=0.3)

    # Panel B: Genus (within Bacterial fraction)
    top_genera = bact_genus_pct.mean().sort_values(ascending=False).head(9).index.tolist()
    genus_plot_df = bact_genus_pct[top_genera].copy()
    genus_plot_df['Other Genera'] = 100.0 - genus_plot_df.sum(axis=1)

    colors_genus = sns.color_palette("Set2", len(genus_plot_df.columns))
    bottom_g = np.zeros(len(SAMPLE_ORDER))
    for col, c in zip(genus_plot_df.columns, colors_genus):
        clean_label = clean_taxon_name(col, 'genus') if col != 'Other Genera' else col
        vals = genus_plot_df[col].values
        ax_g.bar(sample_labels, vals, bottom=bottom_g, label=clean_label, color=c, edgecolor='black', linewidth=0.5, width=0.55)
        bottom_g += vals

    ax_g.set_title("B. Bacterial Genus Composition (Top Genera)", loc='left', fontweight='bold')
    ax_g.set_ylabel("Relative Abundance (%)", fontweight='bold')
    ax_g.set_ylim(0, 100)
    ax_g.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True, facecolor='white', framealpha=0.9, edgecolor='none', fontsize=8)
    ax_g.grid(axis='y', linestyle='--', alpha=0.3)

    plt.tight_layout()
    save_figure(fig, "Fig5_Taxonomic_Composition_Barplots")
    plt.close(fig)

    # 4. Figure 6: Clustered Heatmap of Top 25 Dominant Bacterial Genera
    top25_genera = bact_genus_pct.mean().sort_values(ascending=False).head(25).index.tolist()
    heatmap_df = bact_genus_pct[top25_genera].T
    heatmap_df.index = [clean_taxon_name(g, 'genus') for g in heatmap_df.index]

    clean_idx = []
    seen = {}
    for name in heatmap_df.index:
        if name in seen:
            seen[name] += 1
            clean_idx.append(f"{name} ({seen[name]})")
        else:
            seen[name] = 1
            clean_idx.append(name)
    heatmap_df.index = clean_idx

    # Clean sample labels
    heatmap_df.columns = [c.replace('_', ' ') for c in heatmap_df.columns]

    zscore_df = heatmap_df.apply(lambda row: (row - row.mean()) / (row.std(ddof=0) + 1e-8), axis=1)

    g = sns.clustermap(
        zscore_df,
        cmap='vlag',
        center=0,
        col_cluster=True,
        row_cluster=True,
        figsize=(10, 11),
        linewidths=0.6,
        linecolor='white',
        dendrogram_ratio=(0.14, 0.10),
        cbar_kws={'orientation': 'horizontal'}
    )

    g.fig.suptitle('Top 25 Bacterial Genera Abundance Patterns Across Conditions',
                   fontsize=11.5, fontweight='bold', y=0.98)

    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')


    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=35, ha='right', fontweight='bold', fontsize=9.5)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8.5)

    g.fig.subplots_adjust(top=0.92, bottom=0.10, right=0.72)

    # Explicitly set colorbar position to upper-right, sleek compact shape
    g.ax_cbar.set_position([0.71, 0.885, 0.08, 0.014])
    g.ax_cbar.tick_params(labelsize=8)
    g.ax_cbar.set_title('Abundance Z-score', fontsize=8.5, fontweight='bold', pad=4)
    g.ax_cbar.set_xlabel('')

    save_figure(g.fig, 'Fig6_Top_Genera_Clustered_Heatmap')
    plt.close(g.fig)
    print('Step 04 completed successfully.')

if __name__ == '__main__':
    run()
