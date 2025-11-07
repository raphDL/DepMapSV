PY ?= python
VERBOSE ?= -v
PROGRESS ?= --progress

prep:
	$(PY) prep_depmap_25Q3.py \
		--cnv_segments_wgs OmicsCNSegmentsWGS.csv \
		--expression_tpm OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv \
		--outdir prep_out \
		--cn_assume_diploid \
		$(PROGRESS) $(VERBOSE)
	$(PY) gtf_to_gene_bed.py gencode.v49.basic.annotation.gtf prep_out/gencode.v49.gene.bed

reanalyze:
	$(PY) svbias_reanalysis.py \
		--dependency CRISPRGeneEffect.csv \
		--cnv prep_out/cnv_segments.bed \
		--sv prep_out/sv_from_cnv.bedpe \
			--genes prep_out/genes.depmap.unique.bed \
		--out out \
		--expression prep_out/expression_long.csv \
		--essential CRISPRInferredCommonEssentials.csv \
		--nonessential AchillesNonessentialControls.csv \
		$(PROGRESS) $(VERBOSE)

run: reanalyze

all: prep reanalyze

eval:
	@echo "Evaluation summary (if produced):" && cat out/evaluation_summary.csv || true


