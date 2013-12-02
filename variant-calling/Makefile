export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

PATIENTS = 545 430 460 564 592
DATA = ~/hdall/data
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

all: $(foreach P, $(PATIENTS), $(DATA)/varscan/$P_rem_rel.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.vcf.gz.tbi \
							   $(DATA)/varscan/$P_rem_rel.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.dbsnp.snpeff.vcf.gz.tbi)

%.vcf.gz.tbi: %.vcf
	~/tools/vcftools_0.1.10/bin/vcf-sort $^ | bgzip -c >$*.vcf.gz
	~/tools/tabix-0.2.6/tabix $*.vcf.gz -p vcf	

# -------
# VARSCAN
# -------

$(DATA)/varscan/%_rem_rel.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.vcf: $(DATA)/bam/%_rem.merged.duplicate_marked.realigned.recalibrated.bam $(DATA)/bam/%_rel.merged.duplicate_marked.realigned.recalibrated.bam
	java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar somatic \
		<(~/tools/samtools-0.1.19/samtools view -b -u -q 1 $(word 1,$^) | ~/tools/samtools-0.1.19/samtools mpileup -f ~/generic/data/hg19/ucsc.hg19.fasta -) \
		<(~/tools/samtools-0.1.19/samtools view -b -u -q 1 $(word 2,$^) | ~/tools/samtools-0.1.19/samtools mpileup -f ~/generic/data/hg19/ucsc.hg19.fasta -) \
		$@.part \
		--min-coverage 2 \
		--min-strands2 2 \
		--min-var-freq 0.2 \
		--normal-purity 1 \
		--tumor-purity 0.95 \
		--p-value 1 \
		--somatic-p-value 1 \
		--strand-filter 1 \
		--output-vcf 1 \
		2>&1 | $(LOG)

	bgzip -c $@.part.snp > $@.part.snp.vcf.gz 
	~/tools/tabix-0.2.6/tabix $@.part.snp.vcf.gz -p vcf	
	bgzip -c $@.part.indel > $@.part.indel.vcf.gz 
	~/tools/tabix-0.2.6/tabix $@.part.indel.vcf.gz -p vcf	
	~/tools/vcftools_0.1.10/bin/vcf-concat $@.part.snp.vcf.gz <(~/tools/vcftools_0.1.10/bin/vcf-shuffle-cols -t $@.part.snp.vcf.gz $@.part.indel.vcf.gz) | ~/tools/vcftools_0.1.10/bin/vcf-sort > $@.part
	mv $@.part $@

# relapse and diagnosis samples: more stringent filtering
#$(DATA)/varscan/%_rel.merged.duplicate_marked.realigned.recalibrated.varscan.vcf: $(DATA)/bam/%_rel.merged.duplicate_marked.realigned.recalibrated.bam
#	~/tools/samtools-0.1.19/samtools mpileup -q 40 -f ~/generic/data/hg19/ucsc.hg19.fasta $< \
#		| java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2cns --variants --strand-filter 1 --min-coverage 10 --min-avg-qual 25 --p-value 1 --min-var-freq 0.1 --min-reads2 4 --output-vcf 1 - \
#			> $@.part
#	mv $@.part $@
#	bgzip -c $@ > $@.gz

#$(DATA)/varscan/%_dia.merged.duplicate_marked.realigned.recalibrated.varscan.vcf: $(DATA)/bam/%_dia.merged.duplicate_marked.realigned.recalibrated.bam
#	~/tools/samtools-0.1.19/samtools mpileup -q 40 -f ~/generic/data/hg19/ucsc.hg19.fasta $< \
#		| java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2cns --variants --strand-filter 1 --min-coverage 10 --min-avg-qual 25 --p-value 1 --min-var-freq 0.1 --min-reads2 4 --output-vcf 1 - \
#			> $@.part
#	mv $@.part $@
#	bgzip -c $@ > $@.gz

# remission samples; less stringent filtering
#$(DATA)/varscan/%_rem.merged.duplicate_marked.realigned.recalibrated.varscan.vcf: $(DATA)/bam/%_rem.merged.duplicate_marked.realigned.recalibrated.bam
#	~/tools/samtools-0.1.19/samtools mpileup -q 1 -f ~/generic/data/hg19/ucsc.hg19.fasta $< \
#		| java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2cns --variants --strand-filter 0 --min-coverage 2 --min-avg-qual 15 --p-value 1 --min-var-freq 0.001 --min-reads2 2 --output-vcf 1 - \
#			> $@.part
#	bgzip -c $@ > $@.gz

# -------
# FIND SOMATIC VARIANTS
# -------
 $(DATA)/varscan/%_rel_rem.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.vcf: $(DATA)/varscan/%_rel.merged.duplicate_marked.realigned.recalibrated.varscan.vcf $(DATA)/varscan/%_rem.merged.duplicate_marked.realigned.recalibrated.varscan.vcf ~/hdall/scripts/variant-calling/annotate-remission.pl
	cat $(word 1,$^) \
		| perl ~/hdall/scripts/variant-calling/annotate-remission.pl \
			--patient $* \
			--rem-sample $(word 2,$^).gz \
			> $@.part
	mv $@.part $@ 

 $(DATA)/varscan/%_dia_rem.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.vcf: $(DATA)/varscan/%_dia.merged.duplicate_marked.realigned.recalibrated.varscan.vcf $(DATA)/varscan/%_rem.merged.duplicate_marked.realigned.recalibrated.varscan.vcf ~/hdall/scripts/variant-calling/annotate-remission.pl
	cat $(word 1,$^) \
		| perl ~/hdall/scripts/variant-calling/annotate-remission.pl \
			--patient $* \
			--rem-sample $(word 2,$^).gz \
			> $@.part
	mv $@.part $@ 

# -------
# SNPEFF
# -------

$(DATA)/varscan/%.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.dbsnp.vcf:  $(DATA)/varscan/%.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate \
		-v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf \
		<(cat $< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$@.part | $(LOG))
	mv $@.part $@
 
$(DATA)/varscan/%.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.dbsnp.snpeff.vcf: $(DATA)/varscan/%.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.dbsnp.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.3h; java -Xmx2g -jar snpEff.jar -v -lof hg19 -stats $(DATA)/varscan/$*.snpeff.summary.html $< 2>&1 1>$@.part | $(LOG))
	mv $@.part $@
 
# filter somatic non-silent variants impacting panel genes
$(DATA)/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.snpeff.dbsnp.panel.nonsilent.vcf: $(DATA)/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.snpeff.dbsnp.vcf ~/hdall/scripts/variant-calling/filter-panel-nonsilent.pl 
	cat $< | perl ~/hdall/scripts/variant-calling/filter-panel-nonsilent.pl > $@.part
	mv $@.part $@
	
#-----------	
# FINAL LIST
#-----------	
filtered-variants/%.tsv: $(DATA)/varscan/%.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.dbsnp.snpeff.vcf ~/hdall/scripts/variant-calling/filter-variants.pl
	perl ~/hdall/scripts/variant-calling/filter-variants.pl \
		$< \
		--patient $* \
		--rmsk-file ~/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file ~/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file ~/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file ~/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible ~/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--cosmic-mutation-file ~/generic/data/cosmic/v65/CosmicMutantExport_v65_280513.tsv \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@
	