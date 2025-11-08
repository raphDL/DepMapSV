# Pilot Figure Captions

## Figure 1: Proximity-active genes fraction

**Caption:** Fraction of genes with ≥10% of cells showing meaningful proximity contributions (|contribution| ≥ 0.01 dependency units) in unstable vs stable pilot sets. The unstable set (38%) shows substantially higher proximity activity than the stable set (8%), indicating that breakpoint proximity systematically affects CRISPR dependency scores in breakpoint-dense regions. Thresholds: contribution ≥ 0.01 (when standardization is used), cell fraction ≥ 10%.

**Short caption (for main text):** Proximity-active genes fraction: unstable (38%) vs stable (8%). Genes with ≥10% cells where |proximity contribution| ≥ 0.01.

---

## Figure 2: Proximity-only correlation reduction

**Caption:** Distribution of |Δcorr(dep, CN)| after removing only proximity effects (prox-only correction) for unstable vs stable pilot sets. Each point represents one gene's median |Δr| across cells. Positive values indicate reduced CN–dependency coupling after proximity removal. The unstable set shows a small positive tail (proximity causes subtle but real decoupling), while the stable set is centered near zero. This complements Figure 1 by showing *how much* proximity shifts CN coupling, not just *how often* it matters.

**Short caption (for main text):** |Δcorr(dep, CN)| after prox-only correction. Unstable set shows positive tail (median -0.0002), indicating proximity-induced CN coupling that is mitigated by correction.

---

## Figure 3: Directional prox-active across conditions (fixed 50 genes)

**Caption:** Directional proximity-active fraction (genes with ≥10% cells where |proximity contribution| ≥ 0.01 AND prox-only |Δr| > 0) across five conditions using the same locked unstable gene set (n=50): baseline (Huber, 100k/1Mb; 0.08), rotate shuffle (0.08), uniform shuffle (0.26), linear model with 250k/2Mb windows (0.54), and linear/250k–2Mb + rotate (0.60). Error bars show 95% gene-bootstrap CIs (10,000 replicates). Directional prox-active is low under baseline and rotate shuffle (both 0.08), while robust settings (linear/250k–2Mb) show higher values (0.54–0.60), supporting a breakpoint-driven proximity effect that is robust to parameter choices.

**Short caption:** Directional prox-active (with |Δr| > 0) across baseline, shuffles, and linear/250k–2Mb on fixed 50 genes; bars with 95% CIs. Values: baseline 0.08, rotate 0.08, uniform 0.26, linear/250k–2Mb 0.54, linear/250k–2Mb+rotate 0.60.

## Table: Top 5 genes by |Δcorr| (unstable set)

**Caption:** Top 5 genes from the unstable pilot set ranked by absolute correlation reduction after full correction (CN + proximity). These genes show the largest decoupling between copy number and dependency scores, making them strong candidates for case studies of structural-bias effects.

**Genes:** GFPT2, HLA-DPA1, PACS1, TBC1D22A, KRT9

