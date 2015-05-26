export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
.SECONDEXPANSION:

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PATIENTS = 314 1021247 399 430 446 460 545 564 592 715 786 792 818 842 A B C D E X Y
PATIENTS_XENO = m1963-545-rel m1964-545-rel m1957-715-rel m1977-G-dia m1967-Y-rel

all: circos/allpatients.relapse.circos.png circos/allpatients.diagnosis.circos.png allpatients.filtered-segments.exomeCopy.tsv

counts.RData: /mnt/projects/hdall/data/bam/314_rem.merged.duplicate_marked.realigned.recalibrated.bam
	Rscript /mnt/projects/hdall/scripts/exomeCopy/get-counts.R

counts.bg.RData: counts.RData /mnt/projects/hdall/scripts/exomeCopy/calc-background.R
	Rscript /mnt/projects/hdall/scripts/exomeCopy/calc-background.R

%.combined.pdf %.compiled-segments.tsv: counts.bg.RData /mnt/projects/hdall/scripts/exomeCopy/get-cnv.R
	Rscript /mnt/projects/hdall/scripts/exomeCopy/get-cnv.R --sample $*

allpatients.compiled-segments.exomeCopy.tsv: $(foreach P, $(PATIENTS), $P_rem.compiled-segments.tsv $P_dia.compiled-segments.tsv $P_rel.compiled-segments.tsv) \
											 $(foreach P, $(PATIENTS_XENO), /mnt/projects/p2ry8-crlf2/results/exomeCopy/$P.compiled-segments.tsv) \
											 /mnt/projects/hdall/scripts/exomeCopy/merge-cnvs.R
	Rscript /mnt/projects/hdall/scripts/exomeCopy/merge-cnvs.R
	mv $@.part $@

allpatients.filtered-segments.exomeCopy.tsv: allpatients.compiled-segments.exomeCopy.tsv /mnt/projects/hdall/scripts/exomeCopy/filter-segments.R
	Rscript /mnt/projects/hdall/scripts/exomeCopy/filter-segments.R
	mv $@.part $@
	 
.PHONY: circos
circos: circos/allpatients.diagnosis.circos.png circos/allpatients.relapse.circos.png
 
circos/%.somatic.circos.tsv: allpatients.filtered-segments.exomeCopy.tsv /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/get-circos-somatic.R
	mkdir -p circos
	Rscript /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/get-circos-somatic.R --segments-file $< --output-file $@.part --tumor $* --normal wurscht
	mv $@.part $@ 

circos/allpatients.diagnosis.circos.png: $(foreach P, $(PATIENTS_MATCHED) $(PATIENTS_DIA_ONLY), circos/$PD.somatic.circos.tsv) /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/gene-labels.txt /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/sample-order.dia.txt
	rm -f circos/diagnosis.circos.conf
	perl /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl \
		--template /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf \
		--data-dir /mnt/projects/p2ry8-crlf2/results/exomeCopy/circos \
		--order-file /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/sample-order.dia.txt \
		--file-pattern D.somatic.circos.tsv \
		--gfx-output-file $@.part \
		> circos/diagnosis.circos.conf
	~/tools/circos-0.64/bin/circos --conf circos/diagnosis.circos.conf
	mv $@.part.png $@
	
circos/allpatients.relapse.circos.png: $(foreach P, $(PATIENTS_MATCHED), circos/$PR.somatic.circos.tsv) /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/gene-labels.txt
	rm -f circos/relapse.circos.conf
	perl /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/make-circos-conf.pl \
		--template /mnt/projects/p2ry8-crlf2/scripts/exomeCopy/template.circos.conf \
		--data-dir /mnt/projects/p2ry8-crlf2/results/exomeCopy/circos \
		--file-pattern R.somatic.circos.tsv \
		--gfx-output-file $@.part \
		> circos/relapse.circos.conf
	~/tools/circos-0.64/bin/circos --conf circos/relapse.circos.conf
	mv $@.part.png $@