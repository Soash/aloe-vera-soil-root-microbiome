#!/usr/bin/env python3
"""
Common configuration, publication styling, palettes, and utility loaders
for the Aloe vera 16S Metagenomics Visualization and Analysis Suite.
"""

import os
import sys
import zipfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Color palettes designed for publication (color-blind friendly and distinct)
HEALTH_PALETTE = {
    "Healthy": "#2E7D32",   # Forest Green
    "Disease": "#C62828"    # Crimson Red
}

TYPE_PALETTE = {
    "Root": "#D84315",      # Rust / Root Brown
    "Soil": "#3949AB"       # Deep Indigo / Soil
}

SAMPLE_PALETTE = {
    "Disease_Root": "#C62828",   # Crimson
    "Disease_Soil": "#E65100",   # Deep Orange
    "Healthy_Root": "#2E7D32",   # Dark Green
    "Healthy_Soil": "#1565C0"    # Blue
}

SAMPLE_MARKERS = {
    "Disease_Root": "o",   # Circle
    "Disease_Soil": "s",   # Square
    "Healthy_Root": "^",   # Triangle
    "Healthy_Soil": "D"    # Diamond
}

SAMPLE_ORDER = ["Healthy_Root", "Healthy_Soil", "Disease_Root", "Disease_Soil"]

def set_publication_style():
    """Apply Nature/Cell publication standards to matplotlib rcParams."""
    plt.rcParams.update({
        'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
        'font.family': 'sans-serif',
        'font.size': 9,
        'axes.titlesize': 10,
        'axes.titleweight': 'bold',
        'axes.labelsize': 9,
        'axes.labelweight': 'bold',
        'xtick.labelsize': 8.5,
        'ytick.labelsize': 8.5,
        'legend.fontsize': 8,
        'legend.title_fontsize': 8.5,
        'figure.titlesize': 11,
        'figure.titleweight': 'bold',
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
        'axes.linewidth': 1.0,
        'lines.linewidth': 1.5,
        'patch.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.edgecolor': '#333333',
        'xtick.color': '#333333',
        'ytick.color': '#333333',
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'grid.alpha': 0.3,
        'grid.linestyle': '--'
    })

def save_figure(fig, name, output_dir="results/figures"):
    """Save figure in both high-resolution PNG (300 DPI) and vector PDF formats."""
    os.makedirs(output_dir, exist_ok=True)
    png_path = os.path.join(output_dir, f"{name}.png")
    pdf_path = os.path.join(output_dir, f"{name}.pdf")
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(pdf_path, format='pdf', bbox_inches='tight')
    print(f"  [SAVED FIGURE] {png_path} & {pdf_path}")

def save_table(df, name, output_dir="results/tables", index=True):
    """Save dataframe in both TSV and CSV formats for publications and supplementary data."""
    os.makedirs(output_dir, exist_ok=True)
    tsv_path = os.path.join(output_dir, f"{name}.tsv")
    csv_path = os.path.join(output_dir, f"{name}.csv")
    df.to_csv(tsv_path, sep='\t', index=index)
    df.to_csv(csv_path, sep=',', index=index)
    print(f"  [SAVED TABLE]  {tsv_path} & {csv_path}")

def load_metadata(filepath="metadata.tsv"):
    """Load metadata TSV file."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Metadata file not found: {filepath}")
    df = pd.read_csv(filepath, sep='\t')
    if 'sample-id' in df.columns:
        df = df.set_index('sample-id')
    return df

def extract_from_qzv(qzv_path, inner_pattern):
    """Extract and read a file matching pattern from a QIIME 2 .qzv archive."""
    with zipfile.ZipFile(qzv_path, 'r') as z:
        matches = [f for f in z.namelist() if inner_pattern in f]
        if not matches:
            raise FileNotFoundError(f"No match for '{inner_pattern}' in {qzv_path}")
        target = matches[0]
        return z.open(target)

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
