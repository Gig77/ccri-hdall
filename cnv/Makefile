export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

#PATIENTS = 314 1021247 399 430 446 460 545 592 715 786 792 818 842 A B C D E X Y
PATIENTS = 314 399 430 446 460 545 592 

all: $(foreach P, $(PATIENTS), $P_rem.coverage.tsv $P_dia.coverage.tsv $P_rel.coverage.tsv $P_rem.coverage.mito.tsv $P_dia.coverage.mito.tsv $P_rel.coverage.mito.tsv)

%.coverage.tsv: /data/christian/hdall/data/current/bam/%.merged.duplicate_marked.realigned.recalibrated.bam /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr \
		> $@.part
	mv $@.part $@

%.coverage.mito.tsv: /data/christian/hdall/data/current/bam/%.merged.duplicate_marked.realigned.recalibrated.bam
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-r chrM \
		$< \
	> $@.part
	mv $@.part $@
