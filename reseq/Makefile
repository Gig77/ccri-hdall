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

all: filtered-variants.merged.tsv coverage gene-patient-matrix.reseq.nonrel.tsv filtered-variants.reseq.cosmic.normaf.paper.tsv filtered-variants.reseq.nonrel.paper.tsv gene-patient-matrix.reseq.tsv lolliplot/minaf10/lolliplot_CREBBP_NM_004380_both.svg lolliplot/minaf10/matched/lolliplot_CREBBP_NM_004380_both.svg lolliplot/no-minaf/lolliplot_CREBBP_NM_004380_both.svg

.PHONY: varscan
varscan: $(foreach P, $(PATIENTS_DIA_NON_REL), snpeff/$P_dia.reseq.varscan.dbsnp.snpeff.dbNSFP.vcf)

#-------------
# mutation calling non-relapsing diagnosis samples
#-------------
varscan/%_dia.reseq.varscan.vcf: /mnt/projects/hdall/data/reseq/bam/%_Diagnosis.duplicate_marked.realigned.recalibrated.bam /mnt/projects/hdall/data/reseq/enriched_regions.20padding.bed
	~/tools/samtools-0.1.19/samtools mpileup \
		-q 40 \
		-f /mnt/projects/generic/data/hg19/ucsc.hg19.fasta \
		-l /mnt/projects/hdall/data/reseq/enriched_regions.20padding.bed \
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
snpeff/%.reseq.dbsnp.vcf: /mnt/projects/hdall/data/reseq/somatic_mutations_20bp_padding/reseq_%_calls.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

snpeff/%.reseq.indel.dbsnp.vcf: /mnt/projects/hdall/data/reseq/somatic_indels_20bp_padding/reseq_%_indel.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
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
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar dbnsfp -v /mnt/projects/generic/data/dbNSFP-2.1/dbNSFP2.1.txt $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

#-------------
# lolliplots
#-------------

lolliplot/minaf10/lolliplot_CREBBP_NM_004380_both.svg: ../id-mappings.tsv ../lolliplot/pfam-regions.filtered.tsv filtered-variants.reseq.tsv /mnt/projects/hdall/scripts/lolliplot/lolliplot.pl
	rm -f lolliplot/minaf10/*.svg 
	perl /mnt/projects/hdall/scripts/lolliplot/lolliplot.pl \
		--hugos CREBBP,KRAS,NRAS,PTPN11,FLT3 \
		--filtered-variants filtered-variants.reseq.tsv \
		--min-af 0.1 \
		--output-directory lolliplot/minaf10/ \
		2>&1 | $(LOG)

lolliplot/no-minaf/lolliplot_CREBBP_NM_004380_both.svg: ../id-mappings.tsv ../lolliplot/pfam-regions.filtered.tsv filtered-variants.reseq.tsv /mnt/projects/hdall/scripts/lolliplot/lolliplot.pl
	rm -f lolliplot/no-minaf/*.svg 
	perl /mnt/projects/hdall/scripts/lolliplot/lolliplot.pl \
		--hugos CREBBP,KRAS,NRAS,PTPN11,FLT3 \
		--filtered-variants filtered-variants.reseq.tsv \
		--output-directory lolliplot/no-minaf/ \
		2>&1 | $(LOG)
		
lolliplot/minaf10/matched/lolliplot_CREBBP_NM_004380_both.svg: ../id-mappings.tsv ../lolliplot/pfam-regions.filtered.tsv filtered-variants.reseq.tsv /mnt/projects/hdall/scripts/lolliplot/lolliplot.pl
	mkdir -p lolliplot/minaf10/matched
	rm -f lolliplot/minaf10/matched/*.svg 
	perl /mnt/projects/hdall/scripts/lolliplot/lolliplot.pl \
		--hugos CREBBP,KRAS,NRAS,PTPN11,FLT3 \
		--patients $(subst $(space),$(comma),$(foreach P, $(PATIENTS),$P)) \
		--filtered-variants filtered-variants.reseq.tsv \
		--min-af 0.1 \
		--output-directory lolliplot/minaf10/matched \
		2>&1 | $(LOG)

#-------------
# coverage
#-------------

.PHONY: coverage
coverage: $(foreach P, $(PATIENTS_EXOME), coverage/$P_Diagnosis.coverage.bedtools.txt coverage/$P_Relapse.coverage.bedtools.txt) \
		  $(foreach P, $(PATIENTS_RESEQ_MATCHED), coverage/$P_Diagnosis.coverage.bedtools.txt coverage/$P_Relapse.coverage.bedtools.txt coverage/$P_Remission.coverage.bedtools.txt) \
		  $(foreach P, $(PATIENTS_DIA_ONLY), coverage/$P_Diagnosis.coverage.bedtools.txt coverage/$P_Remission.coverage.bedtools.txt) \
		  $(foreach P, $(PATIENTS_REL_ONLY), coverage/$P_Relapse.coverage.bedtools.txt coverage/$P_Remission.coverage.bedtools.txt) \
		  $(foreach P, $(PATIENTS_DIA_NON_REL), coverage/nonrelapsing/$P_Diagnosis.coverage.bedtools.txt)
		  
coverage/%.coverage.bedtools.txt: /mnt/projects/hdall/data/reseq/bam/%.duplicate_marked.realigned.recalibrated.bam /mnt/projects/hdall/data/reseq/enriched_regions.bed
	samtools view -bq 1 -F 0x400 $< | bedtools coverage -hist -abam - -b /mnt/projects/hdall/data/reseq/enriched_regions.bed | grep ^all > $@.part
	mv $@.part $@

coverage/nonrelapsing/%.coverage.bedtools.txt: /mnt/projects/hdall/data/reseq/bam/%.duplicate_marked.realigned.recalibrated.bam /mnt/projects/hdall/data/reseq/enriched_regions.bed
	samtools view -bq 1 -F 0x400 $< | bedtools coverage -hist -abam - -b /mnt/projects/hdall/data/reseq/enriched_regions.bed | grep ^all > $@.part
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
							 /mnt/projects/hdall/scripts/filter-variants.pl 
	perl  /mnt/projects/hdall/scripts/filter-variants.pl --header >$@.part
	cat filtered_variants/*.indel.filtered.tsv >> $@.part
	cat filtered_variants/*.snp.filtered.tsv >> $@.part
	mv $@.part $@

filtered-variants.reseq.nonrel.tsv: $(foreach P, $(PATIENTS_DIA_NON_REL), filtered_variants/$P_dia.reseq.varscan.filtered.tsv) \
							 		/mnt/projects/hdall/scripts/reseq/filter-variants-nonrel.pl 
	perl  /mnt/projects/hdall/scripts/reseq/filter-variants-nonrel.pl --header >$@.part
	cat filtered_variants/*.varscan.filtered.tsv >> $@.part
	mv $@.part $@

filtered-variants.reseq.nonrel.paper.tsv: filtered-variants.reseq.nonrel.tsv /mnt/projects/hdall/results/panel-genes-paper.tsv
	cat <(head -1 $<) <(grep -wf /mnt/projects/hdall/results/panel-genes-paper.tsv $<) > $@.part
	mv $@.part $@

impacted-genes-list.reseq.nonrel.tsv: filtered-variants.reseq.nonrel.tsv /mnt/projects/hdall/scripts/reseq/impacted-genes-nonrel.pl
	cat $< | perl /mnt/projects/hdall/scripts/reseq/impacted-genes-nonrel.pl \
		2>&1 1>$@.part | tee -a make.log
	mv $@.part $@

gene-patient-matrix.reseq.nonrel.tsv: impacted-genes-list.reseq.nonrel.tsv /mnt/projects/hdall/scripts/reseq/get-gene-patient-matrix-nonrel.pl
	cat $< | perl /mnt/projects/hdall/scripts/reseq/get-gene-patient-matrix-nonrel.pl --mut-details \
		2>&1 1>$@.part | tee -a make.log
	mv $@.part $@

filtered-variants.merged.tsv: filtered-variants.reseq.tsv ../filtered-variants.tsv /mnt/projects/hdall/scripts/reseq/merge_discovery_validation.R
	R --no-save --quiet --slave -f /mnt/projects/hdall/scripts/reseq/merge_discovery_validation.R --args $@.part \
		2>&1 | $(LOG)
	mv $@.part $@
	
filtered-variants.reseq.cosmic.tsv: filtered-variants.reseq.tsv /mnt/projects/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv /mnt/projects/hdall/scripts/annotate-cosmic.pl
	cat $(word 1,$^) | perl /mnt/projects/hdall/scripts/annotate-cosmic.pl \
		--cosmic-mutation-file $(word 2,$^) \
		--only-confirmed \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@ 

filtered-variants.reseq.cosmic.normaf.tsv: filtered-variants.reseq.cosmic.tsv ../cnv/hdall.cnv.tsv /mnt/projects/hdall/scripts/normalize-af.pl
	cat $< | perl /mnt/projects/hdall/scripts/normalize-af.pl \
		--cnv-file ../cnv/hdall.cnv.tsv \
		2>&1 1>$@.part | tee -a make.log
	mv $@.part $@ 

filtered-variants.reseq.cosmic.normaf.paper.tsv: filtered-variants.reseq.cosmic.normaf.tsv /mnt/projects/hdall/results/panel-genes-paper.tsv
	cat <(head -1 $<) <(grep -wf /mnt/projects/hdall/results/panel-genes-paper.tsv $<) > $@.part
	mv $@.part $@
	
filtered_variants/%.reseq.snp.filtered.tsv: snpeff/%.reseq.dbsnp.snpeff.dbNSFP.vcf ../curated-recected-variants.tsv /mnt/projects/hdall/scripts/filter-variants.pl ../remission-variants.tsv.gz.tbi
	perl /mnt/projects/hdall/scripts/filter-variants.pl \
		--sample $* \
		--vcf-in $< \
		--variant-type snp \
		--vcf-out filtered_variants/$*.reseq.dbsnp.snpeff.dbNSFP.filtered.vcf \
		--rmsk-file /mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--rejected-variants-file ../curated-recected-variants.tsv \
		--remission-variants-file ../remission-variants.tsv.gz \
		--evs-file /mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		>$@.part | $(LOG)
	mv $@.part $@

filtered_variants/%.reseq.indel.filtered.tsv: snpeff/%.reseq.indel.dbsnp.snpeff.dbNSFP.vcf ../curated-recected-variants.tsv /mnt/projects/hdall/scripts/filter-variants.pl ../remission-variants.tsv.gz.tbi
	perl /mnt/projects/hdall/scripts/filter-variants.pl \
		--sample $* \
		--vcf-in $< \
		--variant-type indel \
		--vcf-out filtered_variants/$*.reseq.indel.dbsnp.snpeff.dbNSFP.filtered.vcf \
		--rmsk-file /mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--rejected-variants-file ../curated-recected-variants.tsv \
		--remission-variants-file ../remission-variants.tsv.gz \
		--evs-file /mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		>$@.part | $(LOG)
	mv $@.part $@	

filtered_variants/%_dia.reseq.varscan.filtered.tsv: snpeff/%_dia.reseq.varscan.dbsnp.snpeff.dbNSFP.vcf /mnt/projects/hdall/scripts/reseq/filter-variants-nonrel.pl
	perl /mnt/projects/hdall/scripts/reseq/filter-variants-nonrel.pl \
		$< \
		--patient $* \
		--rmsk-file /mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--remission-variants-file /mnt/projects/hdall/results/remission-variants.tsv.gz \
		--cosmic-mutation-file /mnt/projects/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv \
		--evs-file /mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@

#-------------
# PICARD hs_metrics
#-------------

.PHONY: picard
picard: $(foreach P, $(PATIENTS_EXOME), picard/$P_Diagnosis.hs_metrics picard/$P_Relapse.hs_metrics) \
		  $(foreach P, $(PATIENTS_RESEQ_MATCHED), picard/$P_Diagnosis.hs_metrics picard/$P_Relapse.hs_metrics picard/$P_Remission.hs_metrics) \
		  $(foreach P, $(PATIENTS_DIA_ONLY), picard/$P_Diagnosis.hs_metrics picard/$P_Remission.hs_metrics) \
		  $(foreach P, $(PATIENTS_REL_ONLY), picard/$P_Relapse.hs_metrics picard/$P_Remission.hs_metrics) \
		  $(foreach P, $(PATIENTS_DIA_NON_REL), picard/$P_Diagnosis.hs_metrics)

picard/%.hs_metrics: /mnt/projects/hdall/data/reseq/bam/%.duplicate_marked.realigned.recalibrated.bam /mnt/projects/hdall/data/reseq/enriched_regions.noheader.bed
	mkdir -p picard
	~/tools/samtools-0.1.19/samtools view -H $< 2>&1 1> picard/$*.intervallist-picard.bed | $(LOG)
	gawk 'BEGIN { OFS="\t"} {print $$1,$$2,$$3,"+",$$4}' /mnt/projects/hdall/data/reseq/enriched_regions.noheader.bed >> picard/$*.intervallist-picard.bed
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar ~/tools/picard-tools-1.114/CalculateHsMetrics.jar \
		BAIT_INTERVALS=picard/$*.intervallist-picard.bed \
		TARGET_INTERVALS=picard/$*.intervallist-picard.bed \
		INPUT=$< \
		OUTPUT=$@.part \
		REFERENCE_SEQUENCE=/mnt/projects/generic/data/hg19/ucsc.hg19.fasta \
		PER_TARGET_COVERAGE=picard/$*.hs_metrics.per_target_coverage.part \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 | $(LOG)
	mv picard/$*.hs_metrics.per_target_coverage.part picard/$*.hs_metrics.per_target_coverage
	rm picard/$*.intervallist-picard.bed
	mv $@.part $@
	
#--------
# other
#--------

impacted-genes-list.reseq.tsv: filtered-variants.reseq.tsv /mnt/projects/hdall/scripts/impacted-genes.pl
	cat $(word 1,$^) | perl /mnt/projects/hdall/scripts/impacted-genes.pl > $@.part
	mv $@.part $@

gene-patient-matrix.reseq.tsv: impacted-genes-list.reseq.tsv /mnt/projects/hdall/scripts/get-gene-patient-matrix.pl
	cat $(word 1,$^) | perl /mnt/projects/hdall/scripts/get-gene-patient-matrix.pl --mut-details >$@.part
	mv $@.part $@

/mnt/projects/generic/data/dbNSFP-2.1/dbNSFP2.1.txt: 
	#curl http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFPv2.1.zip -o /mnt/projects/generic/data/dbNSFP-2.1/dbNSFPv2.1.zip
	#unzip /mnt/projects/generic/data/dbNSFP-2.1/dbNSFPv2.1.zip
	#(head -n 1 dbNSFP2.1_variant.chr1 ; cat dbNSFP2.1_variant.chr* | grep -v "^#" ) > dbNSFP2.1.txt
	#rm dbNSFP2.1_variant.chr* dbNSFPv2.1.zip

#/mnt/projects/hdall/data/reseq/bam/331test_Diagnosis.duplicate_marked.realigned.recalibrated.bam: /mnt/projects/hdall/data/reseq/bam/331_Diagnosis.duplicate_marked.realigned.recalibrated.bam
#	 ~/tools/samtools-0.1.19/samtools view $< -bL /mnt/projects/hdall/data/reseq/enriched_regions.20padding.bed > $@
