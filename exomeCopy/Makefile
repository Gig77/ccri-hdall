export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
.SECONDEXPANSION:

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PATIENTS = 314 1021247 399 430 446 460 545 564 592 715 786 792 818 842 A B C D E X Y
PATIENTS_XENO = m1963-545-rel m1964-545-rel m1957-715-rel m1977-G-dia m1967-Y-rel

all: allpatients.filtered-segments.exomeCopy.tsv

counts.RData: /mnt/projects/hdall/data/bam/314_rem.merged.duplicate_marked.realigned.recalibrated.bam
	Rscript /mnt/projects/hdall/scripts/exomeCopy/get-counts.R

counts.bg.RData: 
	Rscript /mnt/projects/hdall/scripts/exomeCopy/calc-background.R

%.combined.pdf %.compiled-segments.tsv: counts.bg.RData /mnt/projects/hdall/scripts/exomeCopy/get-cnv.R
	Rscript /mnt/projects/hdall/scripts/exomeCopy/get-cnv.R --sample $*

allpatients.compiled-segments.exomeCopy.tsv: $(foreach P, $(PATIENTS), $P_rem.compiled-segments.tsv $P_dia.compiled-segments.tsv $P_rel.compiled-segments.tsv) \
											 $(foreach P, $(PATIENTS_XENO), $P.compiled-segments.tsv) \
											 /mnt/projects/hdall/scripts/exomeCopy/merge-cnvs.R
	Rscript /mnt/projects/hdall/scripts/exomeCopy/merge-cnvs.R
	mv $@.part $@

allpatients.filtered-segments.exomeCopy.tsv: allpatients.compiled-segments.exomeCopy.tsv /mnt/projects/hdall/scripts/exomeCopy/filter-segments.R
	Rscript /mnt/projects/hdall/scripts/exomeCopy/filter-segments.R
	mv $@.part $@
