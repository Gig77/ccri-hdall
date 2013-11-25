export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

PATIENTS = 545

all: $(foreach S, $(PATIENTS), ~/hdall/data/varscan/$S_rel_rem.merged.duplicate_marked.realigned.recalibrated.varscan.vcf ~/hdall/data/varscan/$S_dia_rem.merged.duplicate_marked.realigned.recalibrated.varscan.vcf)

# download dbSNP data from UCSC
~/generic/data/hg19/hg19.snp137Common.txt: 
	#mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select * from snp137Common' > ~/generic/data/hg19/hg19.snp137Common.txt
	#bgzip -c ~/generic/data/hg19/hg19.snp137Common.txt > ~/generic/data/hg19/hg19.snp137Common.txt.gz
	curl -o ~/generic/data/hg19/hg19.snp137Common.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp137Common.txt.gz
	gunzip ~/generic/data/hg19/hg19.snp137Common.txt.gz
	bgzip -c ~/generic/data/hg19/hg19.snp137Common.txt > ~/generic/data/hg19/hg19.snp137Common.txt.gz
	tabix ~/generic/data/hg19/hg19.snp137Common.txt.gz -s 2 -b 3 -e 4

# relapse and diagnosis samples: more stringent filtering
~/hdall/data/varscan/%.merged.duplicate_marked.realigned.recalibrated.varscan.vcf: ~/hdall/data/bam/%.merged.duplicate_marked.realigned.recalibrated.bam
	~/tools/samtools-0.1.19/samtools mpileup -q 40 -f ~/generic/data/hg19/ucsc.hg19.fasta $< \
		| java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2cns --variants --strand-filter 1 --min-coverage 10 --min-avg-qual 25 --p-value 1 --min-var-freq 0.1 --min-reads2 4 --output-vcf 1 - \
			> $@.part
	mv $@.part $@
	bgzip -c $@ > $@.gz
	tabix -p vcf $@.gz

# remission samples; less stringent filtering
~/hdall/data/varscan/%_rem.merged.duplicate_marked.realigned.recalibrated.varscan.vcf: ~/hdall/data/bam/%_rem.merged.duplicate_marked.realigned.recalibrated.bam
	~/tools/samtools-0.1.19/samtools mpileup -q 1 -f ~/generic/data/hg19/ucsc.hg19.fasta $< \
		| java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2cns --variants --strand-filter 0 --min-coverage 2 --min-avg-qual 15 --p-value 1 --min-var-freq 0.001 --min-reads2 2 --output-vcf 1 - \
			> $@.part
	mv $@.part $@
	bgzip -c $@ > $@.gz
	tabix -p vcf $@.gz

# find somatic variants
 ~/hdall/data/varscan/%_rel_rem.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.vcf: ~/hdall/data/varscan/%_rel.merged.duplicate_marked.realigned.recalibrated.varscan.vcf ~/hdall/data/varscan/%_rem.merged.duplicate_marked.realigned.recalibrated.varscan.vcf.gz ~/hdall/scripts/variant-calling/annotate-remission.pl
	cat $(word 1,$^) \
		| perl ~/hdall/scripts/variant-calling/annotate-remission.pl \
			--patient $* \
			--rem-sample $(word 2,$^).gz \
			> $@.part
	mv $@.part $@ 

 ~/hdall/data/varscan/%_dia_rem.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.vcf: ~/hdall/data/varscan/%_dia.merged.duplicate_marked.realigned.recalibrated.varscan.vcf ~/hdall/data/varscan/%_rem.merged.duplicate_marked.realigned.recalibrated.varscan.vcf.gz ~/hdall/scripts/variant-calling/annotate-remission.pl
	cat $(word 1,$^) \
		| perl ~/hdall/scripts/variant-calling/annotate-remission.pl \
			--patient $* \
			--rem-sample $(word 2,$^).gz \
			> $@.part
	mv $@.part $@ 

# annotate mutation impact on genes using snpEff

~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.snpeff.vcf: ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.vcf
	(cd ~/tools/snpEff-3.3h; java -Xmx2g -jar ~/tools/snpEff-3.3h/snpEff.jar -v -lof hg19 $< > $@.part)
	mv $@.part $@

# annotate dbSnp variants using SnpSift

~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf:
	#cd ~/tools/snpEff-3.3h
	#wget -O dbSnp.vcf.gz ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_no_known_medical_impact_20130930.vcf.gz
	#gunzip common_no_known_medical_impact_20130930.vcf.gz
	#cat common_no_known_medical_impact_20130930.vcf | perl -ne 's/^(\d+|X|Y)\t/chr$1\t/; s/^MT\t/chrM\t/; print $_;' > common_no_known_medical_impact_20130930.chr.vcf
	#rm common_no_known_medical_impact_20130930.vcf
  
~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.snpeff.dbsnp.vcf: ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.snpeff.vcf ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf
	(cd ~/tools/snpEff-3.3h; java -jar SnpSift.jar annotate -v ~/tools/snpEff-3.3h/common_no_known_medical_impact_20130930.chr.vcf $< > $@.part) 
	mv $@.part $@

# filter somatic non-silent variants impacting panel genes
~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.snpeff.dbsnp.panel.nonsilent.vcf: ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.somatic.snpeff.dbsnp.vcf ~/hdall/scripts/variant-calling/filter-panel-nonsilent.pl 
	cat $< | perl ~/hdall/scripts/variant-calling/filter-panel-nonsilent.pl > $@.part
	mv $@.part $@
	
