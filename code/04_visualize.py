#!/usr/bin/env python3
"""Generate integration visualization figures."""

import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse, os


def plot_integration(adata, output_dir="../results"):
    os.makedirs(output_dir, exist_ok=True)

    # ==================== 4-panel overview ====================
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    sc.pl.umap(adata, color="study", ax=axes[0,0], show=False, title="Study",
               legend_loc="right margin", legend_fontsize=7)
    sc.pl.umap(adata, color="platform", ax=axes[0,1], show=False, title="Platform")
    sc.pl.umap(adata, color="cell_type_transferred", ax=axes[1,0], show=False,
               title="Transferred Cell Type", legend_loc="right margin", legend_fontsize=6)
    sc.pl.umap(adata, color="leiden", ax=axes[1,1], show=False,
               title="Leiden Clusters", legend_loc="on data", legend_fontsize=7)
    fig.suptitle(f"DRG Atlas scVI Integration ({adata.n_obs:,} cells)", fontsize=16, y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "integration_overview.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ==================== Label transfer ====================
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    sc.pl.umap(adata, color="cell_type_transferred", ax=axes[0], show=False,
               title="Transferred Labels", legend_loc="right margin", legend_fontsize=7)
    if "transfer_confidence" in adata.obs.columns:
        sc.pl.umap(adata, color="transfer_confidence", ax=axes[1], show=False,
                   title="Confidence", cmap="RdYlGn", vmin=0.3, vmax=1.0)
    else:
        axes[1].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "label_transfer.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ==================== Individual UMAPs ====================
    for col in ["study", "platform", "cell_type_transferred", "leiden"]:
        if col not in adata.obs.columns:
            continue
        fig, ax = plt.subplots(figsize=(8, 7))
        sc.pl.umap(adata, color=col, ax=ax, show=False, title=col)
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, f"umap_{col}.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    # ==================== Organ UMAP (if user data has organ metadata) ====================
    if "organ" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(8, 7))
        sc.pl.umap(adata, color="organ", ax=ax, show=False, title="Target Organ (Retrograde)")
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, "umap_organ.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    # ==================== Cell type composition by study ====================
    if "cell_type_transferred" in adata.obs.columns:
        ct_study = pd.crosstab(adata.obs["cell_type_transferred"], adata.obs["study"])
        ct_norm = ct_study.div(ct_study.sum(axis=1), axis=0)

        fig, ax = plt.subplots(figsize=(10, 8))
        ct_norm.plot(kind="barh", stacked=True, ax=ax)
        ax.set_xlabel("Proportion of cells")
        ax.set_title("Cell Type Composition by Study")
        ax.legend(title="Study", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, "celltype_composition_by_study.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)

    # ==================== Transfer confidence distribution ====================
    if "transfer_confidence" in adata.obs.columns:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Histogram
        axes[0].hist(adata.obs["transfer_confidence"], bins=50,
                     color="steelblue", edgecolor="black", alpha=0.7)
        axes[0].axvline(x=0.5, color="red", linestyle="--", label="0.5 threshold")
        axes[0].set_xlabel("Transfer Confidence")
        axes[0].set_ylabel("Number of cells")
        axes[0].set_title("Label Transfer Confidence Distribution")
        axes[0].legend()

        # Per-study boxplot
        studies = adata.obs["study"].unique()
        conf_data = [adata.obs.loc[adata.obs["study"] == s, "transfer_confidence"].values
                     for s in studies]
        bp = axes[1].boxplot(conf_data, labels=studies, vert=True, patch_artist=True)
        for patch in bp["boxes"]:
            patch.set_facecolor("lightblue")
        axes[1].set_ylabel("Transfer Confidence")
        axes[1].set_title("Confidence by Study")
        axes[1].tick_params(axis="x", rotation=45)

        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, "transfer_confidence_distribution.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)

    print(f"Figures saved to {output_dir}/")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad", help="Path to integrated h5ad")
    parser.add_argument("--output-dir", default="../results")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    plot_integration(adata, args.output_dir)


if __name__ == "__main__":
    main()
