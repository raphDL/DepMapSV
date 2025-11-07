#!/usr/bin/env python3
import os, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from audit_directional import compute_directional_flags

BASE = "out_pilot_baseline_fixed/unstable"
COND_PATHS = {
    "baseline_fixed": "out_pilot_baseline_fixed/unstable",
    "rotate":         "out_pilot_shuffle_rotate/unstable",
    "uniform":        "out_pilot_shuffle_uniform/unstable",
    "linear250k2m":   "out_pilot_lin_250k_2m_fixed/unstable",
    "linear250k1p5m": "out_pilot_lin_250k_1p5m_fixed/unstable",
}

# Optionally include rotate variants if present
opt_2m_rot = "out_pilot_lin_250k_2m_rotate/unstable"
opt_1p5m_rot = "out_pilot_lin_250k_1p5m_rotate/unstable"
if os.path.exists(os.path.join(opt_2m_rot, "design_matrix.csv")):
    COND_PATHS["linear250k2m_rotate"] = opt_2m_rot
if os.path.exists(os.path.join(opt_1p5m_rot, "design_matrix.csv")):
    COND_PATHS["linear250k1p5m_rotate"] = opt_1p5m_rot

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--active-frac", type=float, default=0.10, help="Active cell fraction threshold (default: 0.10)")
    parser.add_argument("--include-excess", action="store_true", help="Also compute and plot excess vs rotate")
    args = parser.parse_args()

    # Lock gene universe (50)
    genes = pd.read_csv(os.path.join(BASE,"design_matrix.csv"))["gene"].drop_duplicates().tolist()
    assert len(genes) == 50, f"Expected 50 genes, got {len(genes)}"

    rows=[]
    for name, path in COND_PATHS.items():
        if not os.path.exists(os.path.join(path, "design_matrix.csv")):
            continue
        flags = compute_directional_flags(path, genes, active_frac_thresh=args.active_frac)
        mean_dir = float(flags["dir_flag"].mean())
        rows.append((name, mean_dir, 50))

    stats = pd.DataFrame(rows, columns=["condition","directional_mean","n_genes"])

    # Bootstrap over genes
    rng = np.random.default_rng(1)
    ci_low, ci_high = [], []
    for _, r in stats.iterrows():
        flags = compute_directional_flags(COND_PATHS[r["condition"]], genes, active_frac_thresh=args.active_frac)["dir_flag"].to_numpy().astype(int)
        bs = [np.mean(flags[rng.integers(0, len(flags), size=len(flags))]) for _ in range(10000)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        ci_low.append(lo)
        ci_high.append(hi)

    stats["ci_low"] = ci_low
    stats["ci_high"] = ci_high

    # Save stats
    os.makedirs("figs_pilot", exist_ok=True)
    stats[["condition","directional_mean","ci_low","ci_high","n_genes"]].to_csv("figs_pilot/final_directional_stats.csv", index=False)

    # Plot
    order = ["baseline_fixed","rotate","uniform","linear250k1p5m","linear250k2m"]
    if "linear250k2m_rotate" in stats["condition"].values:
        order.append("linear250k2m_rotate")
    if "linear250k1p5m_rotate" in stats["condition"].values:
        order.append("linear250k1p5m_rotate")
    plot_df = stats.set_index("condition").loc[[c for c in order if c in stats["condition"].values]].reset_index()

    x = np.arange(len(plot_df))
    y = plot_df["directional_mean"].to_numpy()
    yerr = np.vstack([y - plot_df["ci_low"].to_numpy(), plot_df["ci_high"].to_numpy() - y])

    plt.figure(figsize=(8,3))
    colors = ["#4c78a8","#f58518","#e45756","#54a24b","#72b7b2","#9d755d","#b279a2"]
    plt.bar(x, y, yerr=yerr, capsize=4, color=colors[:len(plot_df)], alpha=0.85)
    plt.xticks(x, [c.replace("_"," ").replace("linear250k2m","linear 250k–2Mb").replace("linear250k1p5m","linear 250k–1.5Mb") for c in plot_df["condition"]], rotation=45, ha="right")
    plt.ylabel(f"Directional prox-active (fraction of 50 genes, active_frac>={args.active_frac})")
    plt.ylim(0, 1.0)
    plt.tight_layout()
    plt.savefig("figs_pilot/final_directional_bar.png", dpi=150)
    plt.savefig("figs_pilot/final_directional_bar.pdf")
    plt.close()

    # Excess vs rotate
    if args.include_excess:
        if "rotate" not in stats["condition"].values:
            print("Warning: rotate condition not found; skipping excess computation")
        else:
            rotate_mean = stats[stats["condition"]=="rotate"]["directional_mean"].iloc[0]
            excess_rows = []
            for _, r in stats.iterrows():
                if r["condition"] == "rotate":
                    continue
                if r["condition"] not in COND_PATHS or not os.path.exists(os.path.join(COND_PATHS[r["condition"]], "design_matrix.csv")):
                    continue
                cond_flags = compute_directional_flags(COND_PATHS[r["condition"]], genes, active_frac_thresh=args.active_frac)["dir_flag"].to_numpy().astype(int)
                rot_flags = compute_directional_flags(COND_PATHS["rotate"], genes, active_frac_thresh=args.active_frac)["dir_flag"].to_numpy().astype(int)
                excess_mean = float(cond_flags.mean() - rot_flags.mean())
                # Paired bootstrap
                bs_excess = []
                for _ in range(10000):
                    idx = rng.integers(0, len(cond_flags), size=len(cond_flags))
                    bs_excess.append(cond_flags[idx].mean() - rot_flags[idx].mean())
                lo, hi = np.percentile(bs_excess, [2.5, 97.5])
                excess_rows.append((r["condition"], excess_mean, lo, hi, 50))
            excess_stats = pd.DataFrame(excess_rows, columns=["condition","excess_mean","ci_low","ci_high","n_genes"])
            excess_stats.to_csv("figs_pilot/final_directional_excess.csv", index=False)

            # Plot excess
            order_excess = [c for c in ["baseline_fixed","uniform","linear250k1p5m","linear250k2m"] if c in excess_stats["condition"].values]
            if "linear250k2m_rotate" in excess_stats["condition"].values:
                order_excess.append("linear250k2m_rotate")
            if "linear250k1p5m_rotate" in excess_stats["condition"].values:
                order_excess.append("linear250k1p5m_rotate")
            plot_excess = excess_stats.set_index("condition").loc[order_excess].reset_index()

            x = np.arange(len(plot_excess))
            y = plot_excess["excess_mean"].to_numpy()
            yerr = np.vstack([y - plot_excess["ci_low"].to_numpy(), plot_excess["ci_high"].to_numpy() - y])

            plt.figure(figsize=(8,3))
            colors_excess = ["#4c78a8","#e45756","#54a24b","#72b7b2","#9d755d","#b279a2"]
            plt.bar(x, y, yerr=yerr, capsize=4, color=colors_excess[:len(plot_excess)], alpha=0.85)
            plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
            plt.xticks(x, [c.replace("_"," ").replace("linear250k2m","linear 250k–2Mb").replace("linear250k1p5m","linear 250k–1.5Mb") for c in plot_excess["condition"]], rotation=45, ha="right")
            plt.ylabel(f"Excess directional vs rotate (fraction of 50 genes, active_frac>={args.active_frac})")
            plt.tight_layout()
            plt.savefig("figs_pilot/final_directional_excess_bar.png", dpi=150)
            plt.savefig("figs_pilot/final_directional_excess_bar.pdf")
            plt.close()
            print("Wrote figs_pilot/final_directional_excess.csv and figs_pilot/final_directional_excess_bar.(png|pdf)")

    print("Wrote figs_pilot/final_directional_stats.csv and figs_pilot/final_directional_bar.(png|pdf)")

if __name__ == "__main__":
    main()
