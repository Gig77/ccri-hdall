export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'
comma:=,
empty:=
space:= $(empty) $(empty)

PATIENTS=1009302 1019357 1020540 1020583 1021247 1021392 1021631 1022914 1023056 1023392 1024518 1024543 1025108 1025409 1187 314 399 430 446 460 545 592 715 786 792 818 842 A B C D FB11 G HS6 K MB2 X 243 933 944 KA14651 KD20493 KE12025 MJ16441 NH17331 PJ13414 RD17412 RS13466 ST13892 ZA16211 Y
PATIENTS_REL_ONLY=1017005 1021865 1023545 AD15 BL16 BM18 CA18 DM1 FE1 FS1 GD18 GD1 HJ15 HJA15 HL1 KA17 KJ17 KL16 LB17 LM18 MJ1 ML10 MV16 NS18 PC16 RT15 RT16 SJM16 SKR1 SL1 SLM1 ST14 WA1 ZE13
PATIENTS_DIA_ONLY=1004564 1010661 1010781 1019964 1020076 1021087 1023338 1023616 1024589 1026233 1026662 B100 EF7 FB14 G44 HD7

all: filtered-variants.merged.tsv filtered-variants.reseq.cosmic.tsv gene-patient-matrix.reseq.tsv lolliplot/lolliplot_CREBBP_NM_004380_both.svg lolliplot/matched/lolliplot_CREBBP_NM_004380_both.svg

#-----------	
# SNPEFF
#-----------	
snpeff/%.reseq.dbsnp.vcf: ~/hdall/data/reseq/somatic_mutations_20bp_padding/reseq_%_calls.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

snpeff/%.reseq.indel.dbsnp.vcf: ~/hdall/data/reseq/somatic_indels_20bp_padding/reseq_%_indel.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

snpeff/%.dbsnp.snpeff.vcf: snpeff/%.dbsnp.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -Xmx2g -jar snpEff.jar -v -lof hg19 -stats $(PWD)/snpeff/$*.snpeff.summary.html $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

snpeff/%.dbsnp.snpeff.dbNSFP.vcf: snpeff/%.dbsnp.snpeff.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar dbnsfp -v ~/generic/data/dbNSFP-2.1/dbNSFP2.1.txt $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

#-------------
# lolliplots
#-------------

lolliplot/lolliplot_CREBBP_NM_004380_both.svg: ../id-mappings.tsv ../lolliplot/pfam-regions.filtered.tsv filtered-variants.reseq.tsv ~/hdall/scripts/lolliplot/lolliplot.pl
	rm -f lolliplot/*.svg 
	perl ~/hdall/scripts/lolliplot/lolliplot.pl \
		--hugos CREBBP,KRAS,NRAS,PTPN11,FLT3 \
		--filtered-variants filtered-variants.reseq.tsv \
		--output-directory lolliplot/ \
		2>&1 | $(LOG)
		
lolliplot/matched/lolliplot_CREBBP_NM_004380_both.svg: ../id-mappings.tsv ../lolliplot/pfam-regions.filtered.tsv filtered-variants.reseq.tsv ~/hdall/scripts/lolliplot/lolliplot.pl
	mkdir -p lolliplot/matched
	rm -f lolliplot/matched/*.svg 
	perl ~/hdall/scripts/lolliplot/lolliplot.pl \
		--hugos CREBBP,KRAS,NRAS,PTPN11,FLT3 \
		--patients $(subst $(space),$(comma),$(foreach P, $(PATIENTS),$P)) \
		--filtered-variants filtered-variants.reseq.tsv \
		--output-directory lolliplot/matched \
		2>&1 | $(LOG)

#-------------
# final lists
#-------------
filtered-variants.reseq.tsv: $(foreach P, $(PATIENTS), filtered_variants/$P_rem_dia.reseq.snp.filtered.tsv) \
							 $(foreach P, $(PATIENTS), filtered_variants/$P_rem_dia.reseq.indel.filtered.tsv) \
							 $(foreach P, $(PATIENTS), filtered_variants/$P_rem_rel.reseq.snp.filtered.tsv) \
							 $(foreach P, $(PATIENTS), filtered_variants/$P_rem_rel.reseq.indel.filtered.tsv) \
							 $(foreach P, $(PATIENTS_DIA_ONLY), filtered_variants/$P_rem_dia.reseq.snp.filtered.tsv) \
							 $(foreach P, $(PATIENTS_DIA_ONLY), filtered_variants/$P_rem_dia.reseq.indel.filtered.tsv) \
							 $(foreach P, $(PATIENTS_REL_ONLY), filtered_variants/$P_rem_rel.reseq.snp.filtered.tsv) \
							 $(foreach P, $(PATIENTS_REL_ONLY), filtered_variants/$P_rem_rel.reseq.indel.filtered.tsv) \
							 filtered_variants/715_rem_rel2.reseq.snp.filtered.tsv \
							 filtered_variants/715_rem_rel2.reseq.indel.filtered.tsv \
							 ~/hdall/scripts/filter-variants.pl 
	perl  ~/hdall/scripts/filter-variants.pl --header >$@.part
	cat filtered_variants/*.filtered.tsv >> $@.part
	mv $@.part $@

filtered-variants.merged.tsv: filtered-variants.reseq.tsv ../filtered-variants.tsv ~/hdall/scripts/reseq/merge_discovery_validation.R
	R --no-save --quiet --slave -f ~/hdall/scripts/reseq/merge_discovery_validation.R --args $@.part \
		2>&1 | $(LOG)
	mv $@.part $@
	
filtered-variants.reseq.cosmic.tsv: filtered-variants.reseq.tsv ~/generic/data/cosmic/v65/CosmicMutantExport_v65_280513.tsv ~/hdall/scripts/annotate-cosmic.pl
	cat $(word 1,$^) | perl ~/hdall/scripts/annotate-cosmic.pl \
		--cosmic-mutation-file $(word 2,$^) \
		--only-confirmed \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@ 

	
filtered_variants/%.reseq.snp.filtered.tsv: snpeff/%.reseq.dbsnp.snpeff.dbNSFP.vcf ../curated-recected-variants.tsv ~/hdall/scripts/filter-variants.pl ../remission-variants.tsv.gz.tbi
	perl ~/hdall/scripts/filter-variants.pl \
		--sample $* \
		--vcf-in $< \
		--variant-type snp \
		--vcf-out filtered_variants/$*.reseq.dbsnp.snpeff.dbNSFP.filtered.vcf \
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

filtered_variants/%.reseq.indel.filtered.tsv: snpeff/%.reseq.indel.dbsnp.snpeff.dbNSFP.vcf ../curated-recected-variants.tsv ~/hdall/scripts/filter-variants.pl ../remission-variants.tsv.gz.tbi
	perl ~/hdall/scripts/filter-variants.pl \
		--sample $* \
		--vcf-in $< \
		--variant-type indel \
		--vcf-out filtered_variants/$*.reseq.indel.dbsnp.snpeff.dbNSFP.filtered.vcf \
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

#--------
# other
#--------

impacted-genes-list.reseq.tsv: filtered-variants.reseq.tsv ~/hdall/scripts/impacted-genes.pl
	cat $(word 1,$^) | perl ~/hdall/scripts/impacted-genes.pl > $@.part
	mv $@.part $@

gene-patient-matrix.reseq.tsv: impacted-genes-list.reseq.tsv ~/hdall/scripts/get-gene-patient-matrix.pl
	cat $(word 1,$^) | perl ~/hdall/scripts/get-gene-patient-matrix.pl --mut-details >$@.part
	mv $@.part $@

~/generic/data/dbNSFP-2.1/dbNSFP2.1.txt: 
	#curl http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFPv2.1.zip -o ~/generic/data/dbNSFP-2.1/dbNSFPv2.1.zip
	#unzip ~/generic/data/dbNSFP-2.1/dbNSFPv2.1.zip
	#(head -n 1 dbNSFP2.1_variant.chr1 ; cat dbNSFP2.1_variant.chr* | grep -v "^#" ) > dbNSFP2.1.txt
	#rm dbNSFP2.1_variant.chr* dbNSFPv2.1.zip