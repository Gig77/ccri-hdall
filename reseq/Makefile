export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'
PATIENTS=1009302 1019357 1020540 1020583 1021247 1021392 1021631 1022914 1023056 1023392 1024518 1024543 1025108 1025409 1187 314 399 430 446 460 545 592 715 786 792 818 842 A B C D FB11 G HS6 K MB2 X

all: filtered-variants.merged.tsv

#-----------	
# SNPEFF
#-----------	
snpeff/%.reseq.dbsnp.snpeff.vcf: snpeff/%.reseq.dbsnp.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -Xmx2g -jar snpEff.jar -v -lof hg19 -stats $(PWD)/snpeff/$*.snpeff.summary.html $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

snpeff/%.reseq.dbsnp.vcf: ~/hdall/data/reseq/mutect_enriched_regions/reseq_%_calls.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

filtered-variants.reseq.tsv: $(foreach P, $(PATIENTS), filtered_variants/$P_rem_dia.reseq.snp.filtered.tsv) \
							 $(foreach P, $(PATIENTS), filtered_variants/$P_rem_rel.reseq.snp.filtered.tsv) \
							 ~/hdall/scripts/filter-variants.pl 
	perl  ~/hdall/scripts/filter-variants.pl --header >$@.part
	cat filtered_variants/*.filtered.tsv >> $@.part
	mv $@.part $@

filtered-variants.merged.tsv: filtered-variants.reseq.tsv ../filtered-variants.tsv ~/hdall/scripts/reseq/merge_discovery_validation.R
	R --no-save --quiet --slave -f ~/hdall/scripts/reseq/merge_discovery_validation.R --args $@.part \
		2>&1 | tee -a make.log
	mv $@.part $@
	
filtered_variants/%.reseq.snp.filtered.tsv: snpeff/%.reseq.dbsnp.snpeff.vcf ../curated-recected-variants.tsv ~/hdall/scripts/filter-variants.pl ../remission-variants.tsv.gz.tbi
	perl ~/hdall/scripts/filter-variants.pl \
		--sample $* \
		--vcf-in $< \
		--variant-type snp \
		--vcf-out filtered_variants/$*.reseq.dbsnp.snpeff.filtered.vcf \
		--rmsk-file ~/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file ~/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file ~/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file ~/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible ~/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro ~/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--rejected-variants-file ../curated-recected-variants.tsv \
		--remission-variants-file ../remission-variants.tsv.gz \
		>$@.part | $(LOG)
	mv $@.part $@
