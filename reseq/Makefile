export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'
comma:=,
empty:=
space:= $(empty) $(empty)

PATIENTS=1009302 1019357 1020540 1020583 1021247 1021392 1021631 1022914 1023056 1023392 1024518 1024543 1025108 1025409 1187 314 399 430 446 460 545 592 715 786 792 818 842 A B C D FB11 G HS6 K MB2 X 243 933 944 KA14651 KD20493 KE12025 MJ16441 NH17331 PJ13414 RD17412 RS13466 ST13892 ZA16211 Y
PATIENTS_EXOME=314 1021247 399 430 446 460 545 592 715 786 792 818 842 A B C D X Y
PATIENTS_RESEQ_MATCHED=1009302 1019357 1020540 1020583 1021392 1021631 1022914 1023056 1023392 1024518 1024543 1025108 1025409 1187 FB11 G HS6 K MB2 243 933 944 KA14651 KD20493 KE12025 MJ16441 NH17331 PJ13414 RD17412 RS13466 ST13892 ZA16211
PATIENTS_REL_ONLY=1017005 1021865 1023545 AD15 BL16 BM18 CA18 DM1 FE1 FS1 GD18 GD1 HJ15 HJA15 HL1 KA17 KJ17 KL16 LB17 LM18 MJ1 ML10 MV16 NS18 PC16 RT15 RT16 SJM16 SKR1 SL1 SLM1 ST14 WA1 ZE13
PATIENTS_DIA_ONLY=1004564 1010661 1010781 1019964 1020076 1021087 1023338 1023616 1024589 1026233 1026662 B100 EF7 FB14 G44 HD7
PATIENTS_DIA_NON_REL=331 380 442 350 461 466 529 591 602 619 633 634 642 653 666 672 697 698 700 709 724 762 776 777 779 782 409 NRD_1 73 NRD_2 NRD_3 60 594 687 748 754 646 530 718 681 39 49 45 54 110 111 134 143 147 199 NRD_4

all: filtered-variants.merged.tsv coverage gene-patient-matrix.reseq.nonrel.tsv filtered-variants.reseq.cosmic.normaf.tsv filtered-variants.reseq.nonrel.tsv gene-patient-matrix.reseq.tsv lolliplot/lolliplot_CREBBP_NM_004380_both.svg lolliplot/matched/lolliplot_CREBBP_NM_004380_both.svg

.PHONY: varscan
varscan: $(foreach P, $(PATIENTS_DIA_NON_REL), snpeff/$P_dia.reseq.varscan.dbsnp.snpeff.dbNSFP.vcf)

#-------------
# mutation calling non-relapsing diagnosis samples
#-------------
varscan/%_dia.reseq.varscan.vcf: ~/hdall/data/reseq/bam/%_Diagnosis.duplicate_marked.realigned.recalibrated.bam ~/hdall/data/reseq/enriched_regions.20padding.bed
	~/tools/samtools-0.1.19/samtools mpileup \
		-q 40 \
		-f ~/generic/data/hg19/ucsc.hg19.fasta \
		-l ~/hdall/data/reseq/enriched_regions.20padding.bed \
		$< \
		| java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2cns \
			--variants \
			--strand-filter 1 \
			--min-coverage 10 \
			--min-avg-qual 15 \
			--p-value 1 \
			--min-var-freq 0.1 \
			--min-reads2 4 \
			--output-vcf 1 - \
			2>&1 1>$@.part | $(LOG)
	cat $@.part | perl -ne 's/\tSample1/\t$*_dia/; print $$_;' > $@
	rm $@.part

snpeff/%.reseq.varscan.dbsnp.vcf: varscan/%.reseq.varscan.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $(PWD)/$< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

#-----------	
# SNPEFF
#-----------	
snpeff/%.reseq.dbsnp.vcf: ~/hdall/data/reseq/somatic_mutations_20bp_padding/reseq_%_calls.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

snpeff/%.reseq.indel.dbsnp.vcf: ~/hdall/data/reseq/somatic_indels_20bp_padding/reseq_%_indel.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
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
# coverage
#-------------

.PHONY: coverage
coverage: $(foreach P, $(PATIENTS_EXOME), coverage/$P_Diagnosis.coverage.bedtools.txt coverage/$P_Relapse.coverage.bedtools.txt) \
		  $(foreach P, $(PATIENTS_RESEQ_MATCHED), coverage/$P_Diagnosis.coverage.bedtools.txt coverage/$P_Relapse.coverage.bedtools.txt coverage/$P_Remission.coverage.bedtools.txt) \
		  $(foreach P, $(PATIENTS_RESEQ_DIAONLY), coverage/$P_Diagnosis.coverage.bedtools.txt coverage/$P_Remission.coverage.bedtools.txt) \
		  $(foreach P, $(PATIENTS_RESEQ_RELONLY), coverage/$P_Relapse.coverage.bedtools.txt coverage/$P_Remission.coverage.bedtools.txt)
		  
coverage/%.coverage.bedtools.txt: ~/hdall/data/reseq/bam/%.duplicate_marked.realigned.recalibrated.bam ~/hdall/data/reseq/enriched_regions.bed
	samtools view -bq 1 -F 0x400 $< | bedtools coverage -hist -abam - -b ~/hdall/data/reseq/enriched_regions.bed | grep ^all > $@.part
	mv $@.part $@

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
	cat filtered_variants/*.indel.filtered.tsv >> $@.part
	cat filtered_variants/*.snp.filtered.tsv >> $@.part
	mv $@.part $@

filtered-variants.reseq.nonrel.tsv: $(foreach P, $(PATIENTS_DIA_NON_REL), filtered_variants/$P_dia.reseq.varscan.filtered.tsv) \
							 		~/hdall/scripts/reseq/filter-variants-nonrel.pl 
	perl  ~/hdall/scripts/reseq/filter-variants-nonrel.pl --header >$@.part
	cat filtered_variants/*.varscan.filtered.tsv >> $@.part
	mv $@.part $@

impacted-genes-list.reseq.nonrel.tsv: filtered-variants.reseq.nonrel.tsv ~/hdall/scripts/reseq/impacted-genes-nonrel.pl
	cat $< | perl ~/hdall/scripts/reseq/impacted-genes-nonrel.pl \
		2>&1 1>$@.part | tee -a make.log
	mv $@.part $@

gene-patient-matrix.reseq.nonrel.tsv: impacted-genes-list.reseq.nonrel.tsv ~/hdall/scripts/reseq/get-gene-patient-matrix-nonrel.pl
	cat $< | perl ~/hdall/scripts/reseq/get-gene-patient-matrix-nonrel.pl --mut-details \
		2>&1 1>$@.part | tee -a make.log
	mv $@.part $@

filtered-variants.merged.tsv: filtered-variants.reseq.tsv ../filtered-variants.tsv ~/hdall/scripts/reseq/merge_discovery_validation.R
	R --no-save --quiet --slave -f ~/hdall/scripts/reseq/merge_discovery_validation.R --args $@.part \
		2>&1 | $(LOG)
	mv $@.part $@
	
filtered-variants.reseq.cosmic.tsv: filtered-variants.reseq.tsv ~/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv ~/hdall/scripts/annotate-cosmic.pl
	cat $(word 1,$^) | perl ~/hdall/scripts/annotate-cosmic.pl \
		--cosmic-mutation-file $(word 2,$^) \
		--only-confirmed \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@ 

filtered-variants.reseq.cosmic.normaf.tsv: filtered-variants.reseq.cosmic.tsv ../cnv/hdall.cnv.tsv ~/hdall/scripts/normalize-af.pl
	cat $< | perl ~/hdall/scripts/normalize-af.pl \
		--cnv-file ../cnv/hdall.cnv.tsv \
		2>&1 1>$@.part | tee -a make.log
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
		--evs-file ~/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
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
		--evs-file ~/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		>$@.part | $(LOG)
	mv $@.part $@	

filtered_variants/%_dia.reseq.varscan.filtered.tsv: snpeff/%_dia.reseq.varscan.dbsnp.snpeff.dbNSFP.vcf ~/hdall/scripts/reseq/filter-variants-nonrel.pl
	perl ~/hdall/scripts/reseq/filter-variants-nonrel.pl \
		$< \
		--patient $* \
		--rmsk-file ~/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file ~/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file ~/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file ~/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible ~/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro ~/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--remission-variants-file ~/hdall/results/remission-variants.tsv.gz \
		--cosmic-mutation-file ~/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv \
		--evs-file ~/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		2>&1 1> $@.part | $(LOG)
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
	
#~/hdall/data/reseq/bam/331test_Diagnosis.duplicate_marked.realigned.recalibrated.bam: ~/hdall/data/reseq/bam/331_Diagnosis.duplicate_marked.realigned.recalibrated.bam
#	 ~/tools/samtools-0.1.19/samtools view $< -bL ~/hdall/data/reseq/enriched_regions.20padding.bed > $@
