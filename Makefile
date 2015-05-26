export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash

PATIENTS = 314 1021247 399 430 446 460 545 592 715 786 792 818 842 A B C D E X Y

all: coverage filtered-variants.tsv filtered-vcf filtered-variants.cosmic.tsv filtered-variants.cosmic.normaf.tsv gene-patient-matrix.tsv gene-patient-matrix.af20.tsv gene-patient-matrix.tier1.tsv cnv/gene-patient-matrix.cnv.tsv lolliplot/lolliplot_CREBBP_NM_004380_both.svg ipa/mutated_relapse.tsv stats

filtered-vcf: $(foreach P, $(PATIENTS), filtered_variants/$P_rem_dia.filtered.vcf.gz filtered_variants/$P_rem_rel.filtered.vcf.gz)

stats: stats/variants-per-chrom-and-ploidy.pdf stats/mutations-per-patient-dia-vs-rel.pdf stats/mutation-profile.af10.pdf stats/variant-coverage.pdf stats/gene-length-bias.pdf clonal-analysis/kernel-densities.allpatients.dia.pdf clonal-analysis/allelic-freq-prob.distributions.pdf clonal-analysis/freq-scatter.pdf clonal-analysis/clonal-progression.pdf

nextera:
	cat kamilla/candidate\ genes\ for\ targeted\ sequencing.tsv | perl /mnt/projects/hdall/scripts/get-nextera-exons.pl --density Standard > kamilla/nextera-exons.standard.csv
	cat kamilla/candidate\ genes\ for\ targeted\ sequencing.tsv | perl /mnt/projects/hdall/scripts/get-nextera-exons.pl --density Dense > kamilla/nextera-exons.dense.csv
	cat kamilla/candidate\ genes\ for\ targeted\ sequencing.tsv | perl /mnt/projects/hdall/scripts/get-nextera-exons.pl --density Standard --cds > kamilla/nextera-exons.standard.cds.csv
	cat kamilla/candidate\ genes\ for\ targeted\ sequencing.tsv | perl /mnt/projects/hdall/scripts/get-nextera-exons.pl --density Dense --cds > kamilla/nextera-exons.dense.cds.csv

filtered-variants.tsv:	$(foreach P, $(PATIENTS), filtered_variants/$P_rem_dia.snp.filtered.tsv) \
						$(foreach P, $(PATIENTS), filtered_variants/$P_rem_dia.indel.filtered.tsv) \
						$(foreach P, $(PATIENTS), filtered_variants/$P_rem_rel.snp.filtered.tsv) \
						$(foreach P, $(PATIENTS), filtered_variants/$P_rem_rel.indel.filtered.tsv) \
						/mnt/projects/hdall/scripts/filter-variants.pl 
	perl  /mnt/projects/hdall/scripts/filter-variants.pl --header 2>&1 1>filtered-variants.tsv.part | tee -a make.log
	cat filtered_variants/*.filtered.tsv >> filtered-variants.tsv.part
	mv filtered-variants.tsv.part filtered-variants.tsv
	
filtered_variants/%.snp.filtered.tsv: /mnt/projects/hdall/data/mutect_vcf/%_calls_snpeff_snpsift.dbsnp.vcf curated-recected-variants.tsv /mnt/projects/hdall/scripts/filter-variants.pl remission-variants.tsv.gz.tbi
	perl /mnt/projects/hdall/scripts/filter-variants.pl \
		--sample $* \
		--vcf-in $< \
		--variant-type snp \
		--vcf-out filtered_variants/$*.snp.filtered.vcf \
		--rmsk-file /mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--rejected-variants-file curated-recected-variants.tsv \
		--remission-variants-file remission-variants.tsv.gz \
		--evs-file /mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		2>&1 1>$@.part | grep -v -P '(Leading or trailing space|variant.Format|on which the variant locates)' | tee -a make.log
	mv $@.part $@

.PRECIOUS: /mnt/projects/hdall/data/mutect_vcf/%_calls_snpeff_snpsift.dbsnp.vcf
/mnt/projects/hdall/data/mutect_vcf/%_calls_snpeff_snpsift.dbsnp.vcf: /mnt/projects/hdall/data/mutect_vcf/%_calls_snpeff_snpsift.vcf /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20130930.chr.vcf
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		> $@.part)
	test -s $@.part
	mv $@.part $@

filtered_variants/%.indel.filtered.tsv: /mnt/projects/hdall/data/somatic_indel_vcf/%_snpeff.dbsnp.vcf curated-recected-variants.tsv /mnt/projects/hdall/scripts/filter-variants.pl remission-variants.tsv.gz.tbi
	perl /mnt/projects/hdall/scripts/filter-variants.pl \
		--sample $* \
		--vcf-in $< \
		--variant-type indel \
		--vcf-out filtered_variants/$*.indel.filtered.vcf \
		--rmsk-file /mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--rejected-variants-file curated-recected-variants.tsv \
		--remission-variants-file remission-variants.tsv.gz \
		--evs-file /mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		2>&1 1>$@.part | grep -v -P '(Leading or trailing space|variant.Format|on which the variant locates)' | tee -a make.log
	mv $@.part $@	
	
filtered_variants/%.filtered.vcf.gz: filtered_variants/%.snp.filtered.vcf filtered_variants/%.indel.filtered.vcf
	bgzip -c $(word 1,$^) > $(word 1,$^).gz 
	~/tools/tabix-0.2.6/tabix $(word 1,$^).gz -p vcf	
	bgzip -c $(word 2,$^) > $(word 2,$^).gz 
	~/tools/tabix-0.2.6/tabix $(word 2,$^).gz -p vcf
	~/tools/vcftools_0.1.10/bin/vcf-concat $(word 1,$^).gz <(~/tools/vcftools_0.1.10/bin/vcf-shuffle-cols -t $(word 1,$^).gz $(word 2,$^).gz) | ~/tools/vcftools_0.1.10/bin/vcf-sort | bgzip -c >$@.part
	mv $@.part $@
		

.PRECIOUS: /mnt/projects/hdall/data/somatic_indel_vcf/%_snpeff.dbsnp.vcf
/mnt/projects/hdall/data/somatic_indel_vcf/%_snpeff.dbsnp.vcf: /mnt/projects/hdall/data/somatic_indel_vcf/%_snpeff.vcf /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20130930.chr.vcf
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate -v /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20130930.chr.vcf $< > $@.part) 
	test -s $@.part
	mv $@.part $@

remission-variants.tsv: /mnt/projects/hdall/data/gatk_vcf_all_patients/analysis_ready_vcf.snpEff.vcf /mnt/projects/hdall/scripts/gatk-vcf2tab.pl
	cat /mnt/projects/hdall/data/gatk_vcf_all_patients/analysis_ready_vcf.snpEff.vcf \
		| perl /mnt/projects/hdall/scripts/gatk-vcf2tab.pl 2>&1 1> $@.part | tee -a make.log
	mv $@.part $@

remission-variants.tsv.gz: remission-variants.tsv
	bgzip -c $^ > $@.part
	mv $@.part $@

remission-variants.tsv.gz.tbi: remission-variants.tsv.gz
	~/tools/tabix-0.2.6/tabix $^ -s 2 -b 3 -e 3

filtered-variants.cosmic.tsv: filtered-variants.tsv /mnt/projects/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv /mnt/projects/hdall/scripts/annotate-cosmic.pl
	cat /mnt/projects/hdall/results/filtered-variants.tsv | perl /mnt/projects/hdall/scripts/annotate-cosmic.pl \
		--cosmic-mutation-file /mnt/projects/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv \
		--only-confirmed \
		2>&1 1>$@.part | tee -a make.log
	mv $@.part $@ 

filtered-variants.cosmic.normaf.tsv: filtered-variants.cosmic.tsv cnv/hdall.cnv.tsv /mnt/projects/hdall/scripts/normalize-af.pl
	cat filtered-variants.cosmic.tsv | perl /mnt/projects/hdall/scripts/normalize-af.pl \
		--cnv-file cnv/hdall.cnv.tsv \
		2>&1 1>$@.part | tee -a make.log
	mv $@.part $@ 

#---
# coverage
#---

.PHONY:
coverage: $(foreach P, $(PATIENTS), coverage/$P_dia.coverage.bedtools.txt coverage/$P_rel.coverage.bedtools.txt coverage/$P_rem.coverage.bedtools.txt) 

coverage/%.coverage.bedtools.txt: /mnt/projects/hdall/data/bam/%.merged.duplicate_marked.realigned.recalibrated.bam /mnt/projects/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed.chr
	samtools view -bq 1 -F 0x400 $< | bedtools coverage -hist -abam - -b /mnt/projects/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed.chr | grep ^all > $@.part
	mv $@.part $@
	
# annotation tracks

/mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt.gz.tbi: /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt
	#mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select * from simpleRepeat' > /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt
	bgzip -c /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt > /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt.gz
	~/tools/tabix-0.2.6/tabix /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt.gz -s 2 -b 3 -e 4

/mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz.tbi: /mnt/projects/generic/data/hg19/hg19.rmsk.txt
	#mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select * from rmsk' > /mnt/projects/generic/data/hg19/hg19.rmsk.txt
	bgzip -c /mnt/projects/generic/data/hg19/hg19.rmsk.txt > /mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz
	~/tools/tabix-0.2.6/tabix /mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz -s 6 -b 7 -e 8

/mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt.gz.tbi: /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt
	#mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select * from genomicSuperDups' > /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt
	bgzip -c /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt > /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt.gz
	~/tools/tabix-0.2.6/tabix /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt.gz -s 2 -b 3 -e 4

/mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz.tbi: /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt
	#mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select * from wgEncodeDacMapabilityConsensusExcludable' > /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt
	bgzip -c /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt > /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz
	~/tools/tabix-0.2.6/tabix /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz -s 2 -b 3 -e 4

/mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz.tbi: /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bb
	#curl -o /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bb http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/paired.end.mapping.1000G..pilot.bb
	#curl -o /mnt/projects/hdall/tools/ucsc/bigBedToBed http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/bigBedToBed
	/mnt/projects/hdall/tools/ucsc/bigBedToBed /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bb /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed
	bgzip -c /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed > /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz
	~/tools/tabix-0.2.6/tabix /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz -s 1 -b 2 -e 3

/mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt.gz.tbi: /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt
	#mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select * from ucscRetroAli5' > /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt
	bgzip -c /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt > /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt.gz
	~/tools/tabix-0.2.6/tabix /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt.gz -s 15 -b 17 -e 18


/mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz.tbi: /mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.txt.tar.gz
	#cd /mnt/projects/generic/data/evs
	#curl -O http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.txt.tar.gz
	#rm -f /mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt*
	#tar xvf ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.txt.tar.gz
	#cat *.txt | grep -vP "^#" | perl -ne 's/^([^:]+):(\S+)\s/chr$1\t$2\t/; print $_;' | sort -t '	' -k1,1 -k2,2n > ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.part
	#rm /mnt/projects/generic/data/evs/*.txt
	#mv ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.part ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt
	bgzip -c ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt > ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz
	~/tools/tabix-0.2.6/tabix ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz -s 1 -b 2
	

impacted-genes-list.tsv: filtered-variants.tsv /mnt/projects/hdall/scripts/impacted-genes.pl
	cat filtered-variants.tsv | perl /mnt/projects/hdall/scripts/impacted-genes.pl \
		2>&1 1>impacted-genes-list.tsv.part | tee -a make.log
	mv impacted-genes-list.tsv.part impacted-genes-list.tsv

# TABLE: filtered-variants
impacted-genes-list.af20.tsv: filtered-variants.tsv /mnt/projects/hdall/scripts/impacted-genes.pl
	cat filtered-variants.tsv \
		| perl -ne 'print $$_ if ((split/\t/)[24] >= 0.20)' \
		| perl /mnt/projects/hdall/scripts/impacted-genes.pl \
		2>&1 1>impacted-genes-list.af20.tsv.part | tee -a make.log
	mv impacted-genes-list.af20.tsv.part impacted-genes-list.af20.tsv

gene-patient-matrix.tsv: impacted-genes-list.tsv /mnt/projects/hdall/scripts/get-gene-patient-matrix.pl
	cat impacted-genes-list.tsv | perl /mnt/projects/hdall/scripts/get-gene-patient-matrix.pl --mut-details \
		--patient-ids 446,314,842,E,Y,460,715,C,399,592,A,786,X,B,430,545,D,792,818,1021247 \
		2>&1 1>gene-patient-matrix.tsv.part | tee -a make.log
	mv gene-patient-matrix.tsv.part gene-patient-matrix.tsv

gene-patient-matrix.af20.tsv: impacted-genes-list.af20.tsv /mnt/projects/hdall/scripts/get-gene-patient-matrix.pl
	cat impacted-genes-list.af20.tsv | perl /mnt/projects/hdall/scripts/get-gene-patient-matrix.pl --mut-details \
		--patient-ids 446,314,842,E,Y,460,715,C,399,592,A,786,X,B,430,545,D,792,818,1021247 \
		2>&1 1>gene-patient-matrix.af20.tsv.part | tee -a make.log
	mv gene-patient-matrix.af20.tsv.part gene-patient-matrix.af20.tsv

# TABLE: filtered-variants
impacted-genes-list.tier1.tsv: filtered-variants.tsv /mnt/projects/hdall/scripts/impacted-genes.pl
	cat filtered-variants.tsv \
		| perl -ne 'print $$_ if (/(patient|HIGH|MODERATE)/ and (split/\t/)[24] >= 0.20)' \
		| perl /mnt/projects/hdall/scripts/impacted-genes.pl \
		2>&1 1>impacted-genes-list.tier1.tsv.part | tee -a make.log
	mv impacted-genes-list.tier1.tsv.part impacted-genes-list.tier1.tsv

gene-patient-matrix.tier1.tsv: impacted-genes-list.tier1.tsv /mnt/projects/hdall/scripts/get-gene-patient-matrix.pl
	cat impacted-genes-list.tier1.tsv | perl /mnt/projects/hdall/scripts/get-gene-patient-matrix.pl --mut-details \
		--patient-ids 446,314,842,E,Y,460,715,C,399,592,A,786,X,B,430,545,D,792,818,1021247 \
		2>&1 1>gene-patient-matrix.tier1.tsv.part | tee -a make.log
	mv gene-patient-matrix.tier1.tsv.part gene-patient-matrix.tier1.tsv

cnv/hdall.cnv.tsv: cnv/CNAsallpatients.tsv cnv/hyperdiploid_CytoHDarray.tsv cnv/715_RR_CytoHD_ChAS.tsv /mnt/projects/hdall/scripts/cnv/parse-cnv-data.pl
	cat cnv/CNAsallpatients.tsv | perl /mnt/projects/hdall/scripts/cnv/parse-cnv-data.pl --format andrea --header \
		2>&1 1>cnv/hdall.cnv.tsv.part | tee -a make.log 
	cat cnv/hyperdiploid_CytoHDarray.tsv | perl /mnt/projects/hdall/scripts/cnv/parse-cnv-data.pl --format maria \
		2>&1 1>>cnv/hdall.cnv.tsv.part | tee -a make.log 
	cat cnv/715_RR_CytoHD_ChAS.tsv | perl /mnt/projects/hdall/scripts/cnv/parse-cnv-data.pl --format maria \
		2>&1 1>>cnv/hdall.cnv.tsv.part | tee -a make.log 
	mv cnv/hdall.cnv.tsv.part cnv/hdall.cnv.tsv

cnv/impacted-genes-list.cnv.tsv: cnv/hdall.cnv.tsv /mnt/projects/hdall/scripts/cnv/impacted-genes-cnv.pl
	cat cnv/hdall.cnv.tsv | perl /mnt/projects/hdall/scripts/cnv/impacted-genes-cnv.pl \
		--max-genes 99999 \
		2>&1 1>cnv/impacted-genes-list.cnv.tsv.part | tee -a make.log 
	mv cnv/impacted-genes-list.cnv.tsv.part cnv/impacted-genes-list.cnv.tsv

cnv/gene-patient-matrix.cnv.tsv: cnv/impacted-genes-list.cnv.tsv /mnt/projects/hdall/scripts/get-gene-patient-matrix.cnv.pl
	cat cnv/impacted-genes-list.cnv.tsv | perl /mnt/projects/hdall/scripts/get-gene-patient-matrix.cnv.pl \
		--max-genes 5 \
		2>&1 1>cnv/gene-patient-matrix.cnv.tsv.part | tee -a make.log
	mv cnv/gene-patient-matrix.cnv.tsv.part cnv/gene-patient-matrix.cnv.tsv

clonal-analysis/allelic-freq-prob.distributions.pdf: /mnt/projects/hdall/scripts/clonal-analysis/estimate-af-dist.R
	R --no-save --quiet --slave -f /mnt/projects/hdall/scripts/clonal-analysis/estimate-af-dist.R \
		2>&1 | tee -a make.log

clonal-analysis/freq-scatter.pdf: /mnt/projects/hdall/results/filtered-variants.cosmic.normaf.tsv /mnt/projects/hdall/scripts/clonal-analysis/freq-scatter.R
	Rscript /mnt/projects/hdall/scripts/clonal-analysis/freq-scatter.R 2>&1

clonal-analysis/kernel-densities.allpatients.dia.pdf: /mnt/projects/hdall/results/filtered-variants.cosmic.normaf.tsv /mnt/projects/hdall/scripts/clonal-analysis/varfreq-histograms.R
	Rscript /mnt/projects/hdall/scripts/clonal-analysis/varfreq-histograms.R 2>&1

clonal-analysis/clonal-progression.pdf: /mnt/projects/hdall/results/filtered-variants.cosmic.normaf.tsv /mnt/projects/hdall/scripts/clonal-analysis/clonal-progression.R
	Rscript /mnt/projects/hdall/scripts/clonal-analysis/clonal-progression.R 2>&1

lolliplot/pfam-regions.filtered.tsv: /mnt/projects/generic/data/pfam-27.0/pfamA.txt /mnt/projects/generic/data/pfam-27.0/Pfam-A.regions.tsv.gz /mnt/projects/hdall/scripts/lolliplot/filter-and-annotate-pfam-regions.pl
	perl /mnt/projects/hdall/scripts/lolliplot/filter-and-annotate-pfam-regions.pl > lolliplot/pfam-regions.filtered.tsv

lolliplot/lolliplot_CREBBP_NM_004380_both.svg: id-mappings.tsv lolliplot/pfam-regions.filtered.tsv filtered-variants.cosmic.tsv /mnt/projects/hdall/scripts/lolliplot/lolliplot.pl
	rm -f lolliplot/*.svg 
	perl /mnt/projects/hdall/scripts/lolliplot/lolliplot.pl \
		--hugos CREBBP,KRAS,NRAS,TRRAP,CDC42EP1,SCN5A,CTBS,FRG1,TBP,USP9X,CACNA1B \
		--filtered-variants filtered-variants.cosmic.tsv \
		--output-directory lolliplot/ \
		2>&1 | tee -a make.log

ipa/mutated_relapse.tsv: music/rel-af20/rem_rel.maf music/rel-af20/smg.tsv
	cut -f 1,9 music/rel-af20/rem_rel.maf | grep -P '(Frame_Shift|In_Frame|Missense|Nonsense|Splice_Site)' | cut -f 1 | sort | uniq | grep -f - music/rel-af20/smg.tsv | cut -f 1,9 | sort -k 2g > ipa/mutated_relapse.tsv.tmp
	grep -vP "^(MESP2|TTN|TBP|CSMD3|DNAH5|RYR1|RYR2|RYR3|DNAH1|DNAH8|DNAH9|MUC2|MUC16|MUC12|MUC5B|OR5H2|OR5H2|OR11H4|OR6F1|OR52R1|OR2T12|OR6V1|OR51I2|OR5I1|OR9Q1|OR9Q1|PDZD7)\t" ipa/mutated_relapse.tsv.tmp > ipa/mutated_relapse.tsv.part
	rm ipa/mutated_relapse.tsv.tmp
	mv ipa/mutated_relapse.tsv.part ipa/mutated_relapse.tsv

ipa/mutated_relapse.noKRASpatients.tsv: music/rel-af20/rem_rel.maf music/rel-af20/smg.tsv
	cut -f 1,9,16 music/rel-af20/rem_rel.maf | grep -P '(1021247|818|842|C|786|X|592|D|399|446|314|792)_rel' | grep -P '(Frame_Shift|In_Frame|Missense|Nonsense|Splice_Site)' | cut -f 1 | sort | uniq | grep -f - music/rel-af20/smg.tsv | cut -f 1,9 | sort -k 2g > ipa/mutated_relapse.noKRASpatients.tsv.tmp
	grep -vP "^(MESP2|TTN|TBP|CSMD3|DNAH5|RYR1|RYR2|RYR3|DNAH1|DNAH8|DNAH9|MUC2|MUC16|MUC12|MUC5B|OR5H2|OR5H2|OR11H4|OR6F1|OR52R1|OR2T12|OR6V1|OR51I2|OR5I1|OR9Q1|OR9Q1|PDZD7)\t" ipa/mutated_relapse.noKRASpatients.tsv.tmp > ipa/mutated_relapse.noKRASpatients.tsv.part
	rm ipa/mutated_relapse.noKRASpatients.tsv.tmp
	mv ipa/mutated_relapse.noKRASpatients.tsv.part ipa/mutated_relapse.noKRASpatients.tsv

ipa/genesize.genelist.ipa.tsv: /mnt/projects/hdall/scripts/pathway-analysis/get-size-genesets.pl
	perl /mnt/projects/hdall/scripts/pathway-analysis/get-size-genesets.pl > $@.part
	mv $@.part $@

stats/variants-per-chrom-and-ploidy.pdf: filtered-variants.cosmic.normaf.tsv /mnt/projects/hdall/scripts/stats/variants-per-chrom-and-ploidy.R
	R --no-save --quiet --slave -f /mnt/projects/hdall/scripts/stats/variants-per-chrom-and-ploidy.R \
		2>&1 | tee -a make.log

stats/mutations-per-patient-dia-vs-rel.pdf: filtered-variants.cosmic.normaf.tsv /mnt/projects/hdall/scripts/stats/mutations-per-patient.R
	R --no-save --quiet --slave -f /mnt/projects/hdall/scripts/stats/mutations-per-patient.R \
		2>&1 | tee -a make.log

stats/mutation-profile.af10.pdf: filtered-variants.tsv /mnt/projects/hdall/scripts/stats/mutation-profile.R
	R --no-save --quiet --slave -f /mnt/projects/hdall/scripts/stats/mutation-profile.R \
		2>&1 | tee -a make.log

stats/variant-coverage.pdf: filtered-variants.cosmic.normaf.tsv /mnt/projects/hdall/scripts/stats/variant-coverage.R
	R --no-save --quiet --slave -f /mnt/projects/hdall/scripts/stats/variant-coverage.R \
		2>&1 | tee -a make.log

stats/gene-length-bias.pdf: gene-size.txt gene-patient-matrix.annotated.tsv /mnt/projects/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed.chr
	R --no-save --quiet --slave -f /mnt/projects/hdall/scripts/stats/gene-length-bias.R \
		2>&1 | tee -a make.log
		
gene-size.txt: id-mappings.tsv /mnt/projects/generic/data/hg19/hg19.knownGene.txt /mnt/projects/hdall/scripts/get-gene-size.pl
	perl /mnt/projects/hdall/scripts/get-gene-size.pl > /mnt/projects/hdall/results/gene-size.txt
	
	