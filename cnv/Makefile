export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

PATIENTS = 314 1021247 399 430 446 460 545 592 715 786 792 818 842 A B C D E X Y
#PATIENTS = 314 399 430 446 460 545 592
PATIENTS_RESEQ = 1024543 430 818

all: $(foreach P, $(PATIENTS), $P_rem.coverage.exome.tsv $P_dia.coverage.exome.tsv $P_rel.coverage.exome.tsv $P_rem.coverage.panel.tsv $P_dia.coverage.panel.tsv $P_rel.coverage.panel.tsv) \
	 $(foreach P, $(PATIENTS_RESEQ), ../reseq/cnv/$P_Remission.coverage.reseq.panel.tsv ../reseq/cnv/$P_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/$P_Relapse.coverage.reseq.panel.tsv)

%.coverage.panel.tsv: /data/christian/hdall/data/current/bam/%.merged.duplicate_marked.realigned.recalibrated.bam /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr ~/git/hdall/cnv/get-exon-coverage.pl
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b ~/hdall/data/reseq/enriched_regions.bed \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed ~/hdall/data/reseq/enriched_regions.bed \
		> $@.part
	mv $@.part $@

%.coverage.exome.tsv: /data/christian/hdall/data/current/bam/%.merged.duplicate_marked.realigned.recalibrated.bam /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr ~/git/hdall/cnv/get-exon-coverage.pl
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr \
		> $@.part
	mv $@.part $@

../reseq/cnv/%.coverage.reseq.panel.tsv: ~/hdall/data/reseq/bam/%.duplicate_marked.realigned.recalibrated.bam ~/hdall/data/reseq/enriched_regions.bed ~/git/hdall/cnv/get-exon-coverage.pl
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b ~/hdall/data/reseq/enriched_regions.bed \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed ~/hdall/data/reseq/enriched_regions.bed \
		> $@.part
	mv $@.part $@

%.coverage.mito.tsv: /data/christian/hdall/data/current/bam/%.merged.duplicate_marked.realigned.recalibrated.bam
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-r chrM \
		$< \
	> $@.part
	mv $@.part $@

cov-plot.%.exome.pdf: %_dia.coverage.exome.tsv %_rel.coverage.exome.tsv %_rem.coverage.exome.tsv ~/hdall/scripts/cnv/cov-plot-exome.R
	Rscript ~/hdall/scripts/cnv/cov-plot-exome.R --patient $* --diagnosis $(word 1,$^) --relapse $(word 2,$^) --remission $(word 3,$^) 2>&1
	mv cov-plot.$*.pdf.part $@

cov-plot.%.panel.pdf: %_dia.coverage.panel.tsv %_rel.coverage.panel.tsv %_rem.coverage.panel.tsv ~/hdall/scripts/cnv/cov-plot-panel.R
	Rscript ~/hdall/scripts/cnv/cov-plot-panel.R --patient $* --diagnosis $(word 1,$^) --relapse $(word 2,$^) --remission $(word 3,$^) 2>&1
	mv cov-plot.$*.reseq.pdf.part $@

../reseq/cnv/cov-plot.%.reseq.panel.pdf: ../reseq/cnv/%_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/%_Relapse.coverage.reseq.panel.tsv ../reseq/cnv/%_Remission.coverage.reseq.panel.tsv ~/hdall/scripts/cnv/cov-plot-panel.R
	Rscript ~/hdall/scripts/cnv/cov-plot-panel.R --patient $* --diagnosis $(word 1,$^) --relapse $(word 2,$^) --remission $(word 3,$^) 2>&1
	mv cov-plot.$*.reseq.pdf.part $@
	
