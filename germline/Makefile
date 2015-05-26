export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

PATIENTS_EXOME = 314 399 430 446 460 545 715 786 792 818 842 A B C D E X Y 1021247 592
PATIENTS_RESEQ_MATCHED = 1009302 1019357 1020540 1020583 1021392 1021631 1022914 1023056 1023392 1024518 1024543 1025108 1025409 1187 FB11 G HS6 K MB2 243 933 944 KA14651 KD20493 KE12025 MJ16441 NH17331 PJ13414 RD17412 RS13466 ST13892 ZA16211
PATIENTS_RESEQ_RELONLY = 1017005 1021865 1023545 AD15 BL16 BM18 CA18 DM1 FE1 FS1 GD18 GD1 HJ15 HJA15 HL1 KA17 KJ17 KL16 LB17 LM18 MJ1 ML10 MV16 NS18 PC16 RT15 RT16 SJM16 SKR1 SL1 SLM1 ST14 WA1 ZE13
PATIENTS_RESEQ_DIAONLY = 1004564 1010661 1010781 1019964 1020076 1021087 1023338 1023616 1024589 1026233 1026662 B100 EF7 FB14 G44 HD7
PATIENTS_RESEQ_NONREL = 331 380 442 350 461 466 529 591 602 619 633 634 642 653 666 672 697 698 700 709 724 762 776 777 779 782 409 NRD_1 73 NRD_2 NRD_3 60 594 687 748 754 646 530 718 681 39 49 45 54 110 111 134 143 147 199 NRD_4

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

all: filtered-variants.germline.nonsilent.tsv

#-----------	
# FINAL LIST
#-----------

filtered-variants.germline.nonsilent.tsv: filtered-variants.germline.tsv /mnt/projects/hdall/scripts/query-text.py
	head -1 $< > $@.part
	python /mnt/projects/hdall/scripts/query-text.py --file $< --sql "select * from data where non_silent = 1" >> $@.part
	mv $@.part $@

filtered-variants.reseq.nonrel.tsv: $(foreach P, $(PATIENTS_DIA_NON_REL), filtered_variants/$P_dia.reseq.varscan.filtered.tsv) \
							 		/mnt/projects/hdall/scripts/reseq/filter-variants-nonrel.pl 
	perl  /mnt/projects/hdall/scripts/reseq/filter-variants-nonrel.pl --header >$@.part
	cat filtered_variants/*.varscan.filtered.tsv >> $@.part
	mv $@.part $@

filtered-variants.germline.tsv: $(foreach P, $(PATIENTS_EXOME), filtered_variants/$P.rel.exome.filtered_variants.germline.tsv) \
								$(foreach P, $(PATIENTS_RESEQ_MATCHED), filtered_variants/$P.dia.reseq.filtered_variants.germline.tsv) \
								$(foreach P, $(PATIENTS_RESEQ_DIAONLY), filtered_variants/$P.dia.reseq.filtered_variants.germline.tsv) \
								$(foreach P, $(PATIENTS_RESEQ_RELONLY), filtered_variants/$P.rel.reseq.filtered_variants.germline.tsv) \
								/mnt/projects/hdall/scripts/germline/filter-variants.pl 
	perl  /mnt/projects/hdall/scripts/germline/filter-variants.pl --header >$@.part
	cat filtered_variants/*.germline.tsv >> $@.part
	mv $@.part $@

filtered_variants/%.filtered_variants.germline.tsv: snpeff/%.varscan2.somatic.dbsnp.snpeff.dbNSFP.vcf /mnt/projects/hdall/scripts/germline/filter-variants.pl
	perl /mnt/projects/hdall/scripts/germline/filter-variants.pl \
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
		--clinvar-file /mnt/projects/generic/data/clinvar/clinvar_20140211.vcf.gz \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@	

# -------
# SNPEFF
# -------

snpeff/%.dbsnp.vcf: vcf/%.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $(PWD)/$< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part | $(LOG))
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

# ------
# copy, merge, sort and index varscan VCF files (provided by knolli for exome sequencing data or created myself)
# ------
vcf/%.rel.exome.varscan2.somatic.vcf: /data/knoll/VarScan2/%.varscan2.somatic.snp.vcf /data/knoll/VarScan2/%.varscan2.somatic.indel.vcf
	bgzip -c $(word 1,$^) > vcf/$*.rel.exome.varscan2.somatic.snp.vcf.gz 
	bgzip -c $(word 2,$^) > vcf/$*.rel.exome.varscan2.somatic.indel.vcf.gz 
	~/tools/vcftools_0.1.10/bin/vcf-concat vcf/$*.rel.exome.varscan2.somatic.snp.vcf.gz vcf/$*.rel.exome.varscan2.somatic.indel.vcf.gz \
		| ~/tools/vcftools_0.1.10/bin/vcf-sort > $@
	bgzip -c $@ > $@.gz 
	~/tools/tabix-0.2.6/tabix $@.gz -p vcf
	test -s $@.gz.tbi

vcf/%.reseq.varscan2.somatic.vcf: varscan/%.reseq.varscan2.somatic.snp.vcf varscan/%.reseq.varscan2.somatic.indel.vcf
	bgzip -c $(word 1,$^) > vcf/$*.reseq.varscan2.somatic.snp.vcf.gz 
	bgzip -c $(word 2,$^) > vcf/$*.reseq.varscan2.somatic.indel.vcf.gz 
	~/tools/vcftools_0.1.10/bin/vcf-concat vcf/$*.reseq.varscan2.somatic.snp.vcf.gz vcf/$*.reseq.varscan2.somatic.indel.vcf.gz \
		| ~/tools/vcftools_0.1.10/bin/vcf-sort > $@
	bgzip -c $@ > $@.gz 
	~/tools/tabix-0.2.6/tabix $@.gz -p vcf
	test -s $@.gz.tbi

# -------
# VARSCAN
# -------

varscan/%.dia.reseq.varscan2.somatic.snp.vcf varscan/%.dia.reseq.varscan2.somatic.indel.vcf: /mnt/projects/hdall/data/reseq/bam/%_Remission.duplicate_marked.realigned.recalibrated.bam /mnt/projects/hdall/data/reseq/bam/%_Diagnosis.duplicate_marked.realigned.recalibrated.bam
	(~/tools/samtools-0.1.19/samtools mpileup -q 1 -f /data/christian/generic/data/current/hg19/ucsc.hg19.fasta $^ \
		| awk '{if($$4 != 0 && $$7 != 0) print $$0}' \
		| java -jar /home/STANNANET/christian.knoll/bin/VarScan2/varscan2-3-6.jar somatic \
			--output-snp varscan/$*.dia.reseq.varscan2.somatic.snp.part \
			--output-indel varscan/$*.dia.reseq.varscan2.somatic.indel.part \
			-mpileup 1 \
			--min-coverage 2 \
			--min-strands2 2 \
			--min-var-freq 0.1 \
			--normal-purity 1 \
			--tumor-purity 0.9 \
			--p-value 1 \
			--somatic-p-value 1 \
			--strand-filter 1 \
			--output-vcf 1) \
		2>&1 | $(LOG)
	mv varscan/$*.dia.reseq.varscan2.somatic.snp.part.vcf varscan/$*.dia.reseq.varscan2.somatic.snp.vcf
	mv varscan/$*.dia.reseq.varscan2.somatic.indel.part.vcf varscan/$*.dia.reseq.varscan2.somatic.indel.vcf
	
varscan/%.rel.reseq.varscan2.somatic.snp.vcf varscan/%.rel.reseq.varscan2.somatic.indel.vcf: /mnt/projects/hdall/data/reseq/bam/%_Remission.duplicate_marked.realigned.recalibrated.bam /mnt/projects/hdall/data/reseq/bam/%_Relapse.duplicate_marked.realigned.recalibrated.bam
	(~/tools/samtools-0.1.19/samtools mpileup -q 1 -f /data/christian/generic/data/current/hg19/ucsc.hg19.fasta $^ \
		| awk '{if($$4 != 0 && $$7 != 0) print $$0}' \
		| java -jar /home/STANNANET/christian.knoll/bin/VarScan2/varscan2-3-6.jar somatic \
			--output-snp varscan/$*.rel.reseq.varscan2.somatic.snp.part \
			--output-indel varscan/$*.rel.reseq.varscan2.somatic.indel.part \
			-mpileup 1 \
			--min-coverage 2 \
			--min-strands2 2 \
			--min-var-freq 0.1 \
			--normal-purity 1 \
			--tumor-purity 0.9 \
			--p-value 1 \
			--somatic-p-value 1 \
			--strand-filter 1 \
			--output-vcf 1) \
		2>&1 | $(LOG)
	mv varscan/$*.rel.reseq.varscan2.somatic.snp.part.vcf varscan/$*.rel.reseq.varscan2.somatic.snp.vcf
	mv varscan/$*.rel.reseq.varscan2.somatic.indel.part.vcf varscan/$*.rel.reseq.varscan2.somatic.indel.vcf