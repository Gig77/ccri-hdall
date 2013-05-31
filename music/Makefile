export SHELLOPTS:=errexit:pipefail

all: dia rel ../gene-patient-matrix.annotated.tsv
dia: ucsc-genes.hg19.roi dia/rem_dia.maf dia/gene_mrs dia/smg.tsv dia/sm_pathways.tsv dia/sm_pathways.annotated.tsv
rel: ucsc-genes.hg19.roi rel/rem_rel.maf rel/gene_mrs rel/smg.tsv rel/sm_pathways.tsv rel/sm_pathways.annotated.tsv

#---

ucsc-genes.hg19.roi: ~/hdall/results/id-mappings.tsv ~/hdall/results/ucsc-genes.hg19.bed ~/hdall/scripts/music/make-roi.pl
	~/tools/lh3-sort/sort -k 1,1N -k 2,2n ~/hdall/results/ucsc-genes.hg19.bed > ~/hdall/results/ucsc-genes.hg19.sorted.bed    # sort UCSC bed file
	~/tools/bedtools-2.17.0/bin/bedtools merge -s -nms -n -i ~/hdall/results/ucsc-genes.hg19.sorted.bed | ~/tools/lh3-sort/sort -k 1,1N -k 2,2n > ~/hdall/results/ucsc-genes.hg19.sorted.merged.bed   # merge overlapping exons using bedtools

	# postprocess ROI (map back to gene symbols, re-map some IDs, 1-based start coordinate, cut back coordinates exceeding chromosome length)
	cat ~/hdall/results/ucsc-genes.hg19.sorted.merged.bed | perl ~/hdall/scripts/music/make-roi.pl > ucsc-genes.hg19.roi.part
	mv ucsc-genes.hg19.roi.part ucsc-genes.hg19.roi

	# compress and index it for fast access with tabix
	bgzip -c ucsc-genes.hg19.roi > ucsc-genes.hg19.roi.gz
	tabix ucsc-genes.hg19.roi.gz -s 1 -b 2 -e 3

#---

dia/rem_dia.maf: ~/hdall/results/filtered_variants/314_rem_dia.snp.filtered.vcf ucsc-genes.hg19.roi 

PATIENTS = 314 1021247 399 430 446 460 545 592 715 786 792 818 842 A B C D E X Y

dia/rem_dia.maf: $(foreach P, $(PATIENTS), maf/$P_rem_dia.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>dia/rem_dia.maf.part | tee -a make.log
	cat maf/*rem_dia.maf >> dia/rem_dia.maf.part
	mv dia/rem_dia.maf.part dia/rem_dia.maf

rel/rem_rel.maf: $(foreach P, $(PATIENTS), maf/$P_rem_rel.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>rel/rem_rel.maf.part | tee -a make.log
	cat maf/*rem_rel.maf >> rel/rem_rel.maf.part
	mv rel/rem_rel.maf.part rel/rem_rel.maf

maf/%.maf: ../filtered_variants/%.snp.filtered.vcf ../filtered_variants/%.indel.filtered.vcf ucsc-genes.hg19.roi ~/hdall/results/gene-id-mapping.biomart-0.7.tsv
	mkdir -p maf
	cat ../filtered_variants/$*.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		2>&1 1>$@.part | grep -vP '(Leading or trailing space|variant.Format)' | tee -a make.log
	cat ../filtered_variants/$*.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		2>&1 1>>$@.part | grep -vP '(Leading or trailing space|variant.Format)' | tee -a make.log
	mv $@.part $@
		
#---
	
dia/total_covgs: dia/wig-list.tsv ~/hdall/data/hg19/ucsc.hg19.nochr.fasta
	rm -f dia/total_covgs
	rm -rf dia/roi_covgs
	rm -rf dia/gene_covgs
	genome music bmr calc-wig-covg \
		--wig-list dia/wig-list.tsv \
		--output-dir dia/ \
		--reference-sequence ~/hdall/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

rel/total_covgs: rel/wig-list.tsv ~/hdall/data/hg19/ucsc.hg19.nochr.fasta
	rm -f rel/total_covgs
	rm -rf rel/roi_covgs
	rm -rf rel/gene_covgs
	genome music bmr calc-wig-covg \
		--wig-list rel/wig-list.tsv \
		--output-dir rel/ \
		--reference-sequence ~/hdall/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

dia+rel/total_covgs: dia/total_covgs rel/total_covgs
	rm -f dia+rel/total_covgs
	rm -rf dia+rel/roi_covgs
	rm -rf dia+rel/gene_covgs
	cat dia/total_covgs rel/total_covgs | perl -ne '$$line++; print $$_ if ($$line == 1 or !/^#/);' > dia+rel/total_covgs
	cp -r dia/roi_covgs dia+rel/roi_covgs
	cp rel/roi_covgs/* dia+rel/roi_covgs/
	cp -r dia/gene_covgs dia+rel/gene_covgs
	cp rel/gene_covgs/* dia+rel/gene_covgs/

#---

dia/gene_mrs: dia/total_covgs dia/wig-list.tsv dia/rem_dia.maf ~/hdall/data/hg19/ucsc.hg19.nochr.fasta ucsc-genes.hg19.roi
	rm -f dia/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list dia/wig-list.tsv \
		--maf-file dia/rem_dia.maf \
		--output-dir dia/ \
		--reference-sequence ~/hdall/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

rel/gene_mrs: rel/total_covgs rel/wig-list.tsv rel/rem_rel.maf ~/hdall/data/hg19/ucsc.hg19.nochr.fasta ucsc-genes.hg19.roi
	rm -f rel/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list rel/wig-list.tsv \
		--maf-file rel/rem_rel.maf \
		--output-dir rel/ \
		--reference-sequence ~/hdall/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

dia+rel/gene_mrs: dia+rel/total_covgs dia+rel/wig-list.tsv dia+rel/rem_dia_rel.maf ~/hdall/data/hg19/ucsc.hg19.nochr.fasta ucsc-genes.hg19.roi
	rm -f dia+rel/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list dia+rel/wig-list.tsv \
		--maf-file dia+rel/rem_dia_rel.maf \
		--output-dir dia+rel/ \
		--reference-sequence ~/hdall/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

#--- GENOME MUSIC SMG

dia/smg.tsv: dia/gene_mrs
	genome music smg \
		--gene-mr-file dia/gene_mrs \
		--output-file dia/smg.tsv.part
	mv dia/smg.tsv.part dia/smg.tsv
		
rel/smg.tsv: rel/gene_mrs
	genome music smg \
		--gene-mr-file rel/gene_mrs \
		--output-file rel/smg.tsv.part
	mv rel/smg.tsv.part rel/smg.tsv

pathways-david-complete.tsv: pathways-david.tsv pathways-david-reverse.tsv
	perl ~/hdall/scripts/merge-david-annotations.pl \
		pathways-david.tsv \
		pathways-david-reverse.tsv \
		> pathways-david-complete.tsv.part
	mv pathways-david-complete.tsv.part pathways-david-complete.tsv

#--- GENOME MUSIC PATH-SCAN

dia/sm_pathways.tsv: dia/smg.tsv pathways-david-complete.tsv
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia/gene_covgs/ \
		--maf-file dia/rem_dia.maf \
		--output-file dia/sm_pathways.tsv.part \
		--pathway-file pathways-david-complete.tsv
	mv dia/sm_pathways.tsv.part dia/sm_pathways.tsv

dia/sm_pathways.annotated.tsv: dia/sm_pathways.tsv dia/sm_pathways.tsv_detailed
	perl ~/hdall/scripts/music/annotate-pathscan-result.pl \
		--sm-pathways dia/sm_pathways.tsv \
		--sm-pathways-detail dia/sm_pathways.tsv_detailed \
		> dia/sm_pathways.annotated.tsv.part
	mv dia/sm_pathways.annotated.tsv.part dia/sm_pathways.annotated.tsv

rel/sm_pathways.tsv: rel/smg.tsv pathways-david-complete.tsv
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel/gene_covgs/ \
		--maf-file rel/rem_rel.maf \
		--output-file rel/sm_pathways.tsv.part \
		--pathway-file pathways-david-complete.tsv
	mv rel/sm_pathways.tsv.part rel/sm_pathways.tsv

rel/sm_pathways.annotated.tsv: rel/sm_pathways.tsv rel/sm_pathways.tsv_detailed
	perl ~/hdall/scripts/music/annotate-pathscan-result.pl \
		--sm-pathways rel/sm_pathways.tsv \
		--sm-pathways-detail rel/sm_pathways.tsv_detailed \
		> rel/sm_pathways.annotated.tsv.part
	mv rel/sm_pathways.annotated.tsv.part rel/sm_pathways.annotated.tsv
	
#--- COMBINE ALL RESULTS INTO FINAL LIST

../gene-patient-matrix.annotated.tsv: ../gene-patient-matrix.details.tsv dia/smg.tsv rel/smg.tsv dia/sm_pathways.annotated.tsv rel/sm_pathways.annotated.tsv ~/hdall/scripts/annotate-gene-patient-matrix.pl
	perl ~/hdall/scripts/annotate-gene-patient-matrix.pl \
		--gene-patient-matrix ../gene-patient-matrix.details.tsv \
		--smg-dia dia/smg.tsv \
		--smg-rel rel/smg.tsv \
		--smp-dia dia/sm_pathways.annotated.tsv \
		--smp-rel rel/sm_pathways.annotated.tsv \
		2>&1 1>../gene-patient-matrix.annotated.tsv.part | tee -a make.log
	mv ../gene-patient-matrix.annotated.tsv.part ../gene-patient-matrix.annotated.tsv	