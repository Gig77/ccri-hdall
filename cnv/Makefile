export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

#PATIENTS = 314 1021247 399 430 446 460 545 592 715 786 792 818 842 A B C D E X Y
#PATIENTS = 314 399 430 446 460 545 592 
PATIENTS_RESEQ = 1024543 430 818

all: $(foreach P, $(PATIENTS), $P_rem.coverage.tsv $P_dia.coverage.tsv $P_rel.coverage.tsv $P_rem.coverage.mito.tsv $P_dia.coverage.mito.tsv $P_rel.coverage.mito.tsv cov-plot.$P.pdf) \
	 $(foreach P, $(PATIENTS_RESEQ), $P_Remission.coverage.reseq.tsv $P_Diagnosis.coverage.reseq.tsv $P_Relapse.coverage.reseq.tsv cov-plot.$P.reseq.pdf)

%.coverage.tsv: /data/christian/hdall/data/current/bam/%.merged.duplicate_marked.realigned.recalibrated.bam /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr ~/git/hdall/cnv/get-exon-coverage.pl
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr \
		> $@.part
	mv $@.part $@

%.coverage.reseq.tsv: ~/hdall/data/reseq/bam/%.duplicate_marked.realigned.recalibrated.bam ~/hdall/data/reseq/enriched_regions.bed ~/git/hdall/cnv/get-exon-coverage.pl
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

cov-plot.%.pdf: %_dia.coverage.tsv %_rel.coverage.tsv %_rem.coverage.tsv ~/hdall/scripts/cnv/cov-plot.R
	R --no-save --quiet --slave -f ~/hdall/scripts/cnv/cov-plot.R --args $* $(word 1,$^) $(word 2,$^) $(word 3,$^) 2>&1
	mv $@.part $@

cov-plot.%.reseq.pdf: %_Diagnosis.coverage.reseq.tsv %_Relapse.coverage.reseq.tsv %_Remission.coverage.reseq.tsv ~/hdall/scripts/cnv/cov-plot.reseq.R
	R --no-save --quiet --slave -f ~/hdall/scripts/cnv/cov-plot.reseq.R --args $* $(word 1,$^) $(word 2,$^) $(word 3,$^) 2>&1
	mv $@.part $@
	
