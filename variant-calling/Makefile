export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work

# download dbSNP data from UCSC
~/hdall/data/hg19/hg19.snp137Common.txt: 
	#mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select * from snp137Common' > ~/hdall/data/hg19/hg19.snp137Common.txt
	#bgzip -c ~/hdall/data/hg19/hg19.snp137Common.txt > ~/hdall/data/hg19/hg19.snp137Common.txt.gz
	curl -o ~/hdall/data/hg19/hg19.snp137Common.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp137Common.txt.gz
	gunzip ~/hdall/data/hg19/hg19.snp137Common.txt.gz
	bgzip -c ~/hdall/data/hg19/hg19.snp137Common.txt > ~/hdall/data/hg19/hg19.snp137Common.txt.gz
	tabix ~/hdall/data/hg19/hg19.snp137Common.txt.gz -s 2 -b 3 -e 4

varscan: ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.snp ~/hdall/data/varscan/545_rem.merged.duplicate_marked.realigned.recalibrated.varscan.snp

~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.snp: ~/hdall/data/bam/545_rel.merged.duplicate_marked.realigned.recalibrated.bam
	~/tools/samtools-0.1.19/samtools mpileup \
		-q 1 \
		-f ~/hdall/data/hg19/ucsc.hg19.fasta \
		~/hdall/data/bam/545_rel.merged.duplicate_marked.realigned.recalibrated.bam \
		| tee >(java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2snp - > ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.snp.part) \
		| java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2indel - > ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.indel.part
	mv ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.snp.part ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.snp
	mv ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.indel.part ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.indel

~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.vcf: ~/hdall/data/bam/545_rel.merged.duplicate_marked.realigned.recalibrated.bam
	~/tools/samtools-0.1.19/samtools mpileup \
		-q 1 \
		-f ~/hdall/data/hg19/ucsc.hg19.fasta \
		~/hdall/data/bam/545_rel.merged.duplicate_marked.realigned.recalibrated.bam \
		| java -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2cns \
			--variants \
			--p-value 1 \
			--output-vcf 1 - \
			> ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.vcf.part
	mv ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.vcf.part ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.vcf

~/hdall/results/variant-calling/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.snpeff.vcf:
	java -Xmx2g -jar ~/tools/snpEff-3.3h/snpEff.jar -v -lof -hgvs hg19 ~/hdall/data/varscan/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.vcf > ~/hdall/results/variant-calling/545_rel.merged.duplicate_marked.realigned.recalibrated.varscan.snpeff.vcf 