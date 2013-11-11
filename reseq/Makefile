export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'
PATIENTS=1009302 1019357 1020540 1020583 1021247 1021392 1021631 1022914 1023056 1023392 1024518 1024543 1025108 1025409 1187 314 399 430 446 460 545 592 715 786 792 818 842 A B C D FB11 G HS6 K MB2 X

all: $(foreach P, $(PATIENTS), snpeff/$P_rem_dia.reseq.dbsnp.snpeff.vcf) $(foreach P, $(PATIENTS), snpeff/$P_rem_rel.reseq.dbsnp.snpeff.vcf)

#-----------	
# SNPEFF
#-----------	
snpeff/%.reseq.dbsnp.snpeff.vcf: snpeff/%.reseq.dbsnp.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -Xmx2g -jar snpEff.jar -v -lof hg19 $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv ~/tools/snpEff-3.3h/snpEff_summary.html snpeff/$*.snpEff_summary.html
	mv ~/tools/snpEff-3.3h/snpEff_genes.txt snpeff/$*.snpEff_genes.txt
	mv $@.part $@

snpeff/%.reseq.dbsnp.vcf: ~/hdall/data/reseq/mutect_enriched_regions/reseq_%_calls.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate -v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf $< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@
