#!/bin/bash
# Commands to run focused pilot on top excess candidate genes

cd /Users/raphael/Documents/sandbox/DepMapSV

echo "=========================================="
echo "FOCUSED PILOT: Top Excess Candidate Genes"
echo "=========================================="
echo ""
echo "Selected genes: $(wc -l < figs_v2/pilot_genes.txt) genes"
echo ""

echo "Run focused pilot on top candidates (TRUE):"
echo ""
echo "python svbias_reanalysis.py --pilot --pilot-genes-file figs_v2/pilot_genes.txt \\"
echo "  --dependency CRISPRGeneEffect.csv --cnv prep_out/cnv_segments.bed --sv prep_out/sv_from_cnv.bedpe \\"
echo "  --genes prep_out/genes.depmap.unique.bed --model linear --bp_windows 250000 2000000 \\"
echo "  --add-continuous-proximity --standardize-predictors --out out_focus_true --progress -vv"
echo ""

echo "Run rotate-shuffle control on same genes:"
echo ""
echo "python svbias_reanalysis.py --pilot --pilot-genes-file figs_v2/pilot_genes.txt --shuffle-rotate \\"
echo "  --dependency CRISPRGeneEffect.csv --cnv prep_out/cnv_segments.bed --sv prep_out/sv_from_cnv.bedpe \\"
echo "  --genes prep_out/genes.depmap.unique.bed --model linear --bp_windows 250000 2000000 \\"
echo "  --add-continuous-proximity --standardize-predictors --out out_focus_rotate --progress -vv"
echo ""

echo "=========================================="
echo "SUCCESS/FAILURE CRITERIA:"
echo "=========================================="
echo ""
echo "If TRUE beats rotate on:"
echo "  - Directional prox-active fraction"
echo "  - Median Δ|corr| for selected genes"
echo ""
echo "→ Model-spec issue (richer spec helps real cases)"
echo ""
echo "If rotate still ≥ TRUE:"
echo "→ Shuffled worlds are easier to 'fix' (simpler than reality)"
echo "→ Consider narrow-scope paper or pivot"
echo ""

