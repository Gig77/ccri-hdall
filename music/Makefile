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

dia/rem_dia.maf: ~/hdall/results/filtered_vcf/314_rem_dia.snp.filtered.vcf ucsc-genes.hg19.roi 

	cat ~/hdall/results/filtered_vcf/314_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 314_dia --sample-normal 314_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv --header > dia/rem_dia.maf 2>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/314_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 314_dia --sample-normal 314_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/399_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 399_dia --sample-normal 399_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/399_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 399_dia --sample-normal 399_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/430_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 430_dia --sample-normal 430_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/430_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 430_dia --sample-normal 430_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/446_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 446_dia --sample-normal 446_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/446_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 446_dia --sample-normal 446_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/460_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 460_dia --sample-normal 460_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/460_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 460_dia --sample-normal 460_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/545_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 545_dia --sample-normal 545_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/545_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 545_dia --sample-normal 545_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/592_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 592_dia --sample-normal 592_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/592_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 592_dia --sample-normal 592_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/715_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 715_dia --sample-normal 715_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/715_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 715_dia --sample-normal 715_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/786_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 786_dia --sample-normal 786_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/786_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 786_dia --sample-normal 786_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/792_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 792_dia --sample-normal 792_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/792_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 792_dia --sample-normal 792_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/818_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 818_dia --sample-normal 818_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/818_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 818_dia --sample-normal 818_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/842_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 842_dia --sample-normal 842_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/842_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 842_dia --sample-normal 842_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/1021247_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 1021247_dia --sample-normal 1021247_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/1021247_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 1021247_dia --sample-normal 1021247_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/A_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor A_dia --sample-normal A_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/A_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor A_dia --sample-normal A_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/B_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor B_dia --sample-normal B_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/B_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor B_dia --sample-normal B_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/C_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor C_dia --sample-normal C_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/C_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor C_dia --sample-normal C_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/D_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor D_dia --sample-normal D_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/D_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor D_dia --sample-normal D_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/E_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor E_dia --sample-normal E_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/E_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor E_dia --sample-normal E_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/X_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor X_dia --sample-normal X_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/X_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor X_dia --sample-normal X_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/Y_rem_dia.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor Y_dia --sample-normal Y_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log
	cat ~/hdall/results/filtered_vcf/Y_rem_dia.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor Y_dia --sample-normal Y_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> dia/rem_dia.maf 2>>dia/rem_dia.maf.log

rel/rem_rel.maf: ~/hdall/results/filtered_vcf/314_rem_rel.snp.filtered.vcf ucsc-genes.hg19.roi

	cat ~/hdall/results/filtered_vcf/314_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 314_rel --sample-normal 314_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv --header > rel/rem_rel.maf 2>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/314_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 314_rel --sample-normal 314_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/399_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 399_rel --sample-normal 399_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/399_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 399_rel --sample-normal 399_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/430_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 430_rel --sample-normal 430_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/430_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 430_rel --sample-normal 430_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/446_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 446_rel --sample-normal 446_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/446_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 446_rel --sample-normal 446_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/460_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 460_rel --sample-normal 460_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/460_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 460_rel --sample-normal 460_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/545_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 545_rel --sample-normal 545_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/545_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 545_rel --sample-normal 545_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/592_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 592_rel --sample-normal 592_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/592_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 592_rel --sample-normal 592_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/715_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 715_rel --sample-normal 715_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/715_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 715_rel --sample-normal 715_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/786_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 786_rel --sample-normal 786_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/786_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 786_rel --sample-normal 786_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/792_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 792_rel --sample-normal 792_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/792_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 792_rel --sample-normal 792_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/818_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 818_rel --sample-normal 818_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/818_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 818_rel --sample-normal 818_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/842_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 842_rel --sample-normal 842_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/842_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 842_rel --sample-normal 842_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/1021247_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 1021247_rel --sample-normal 1021247_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/1021247_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor 1021247_rel --sample-normal 1021247_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/A_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor A_rel --sample-normal A_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/A_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor A_rel --sample-normal A_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/B_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor B_rel --sample-normal B_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/B_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor B_rel --sample-normal B_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/C_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor C_rel --sample-normal C_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/C_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor C_rel --sample-normal C_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/D_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor D_rel --sample-normal D_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/D_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor D_rel --sample-normal D_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/E_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor E_rel --sample-normal E_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/E_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor E_rel --sample-normal E_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/X_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor X_rel --sample-normal X_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/X_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor X_rel --sample-normal X_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/Y_rem_rel.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor Y_rel --sample-normal Y_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log
	cat ~/hdall/results/filtered_vcf/Y_rem_rel.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl --sample-tumor Y_rel --sample-normal Y_rem --music-roi ucsc-genes.hg19.roi.gz --mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv >> rel/rem_rel.maf 2>>rel/rem_rel.maf.log

dia+rel/rem_dia_rel.maf: dia/rem_dia.maf rel/rem_rel.maf
	# combine MAFs, removing duplicates
	cat dia/rem_dia.maf rel/rem_rel.maf | perl -ne '@f = split/\t/; print "$_" if (!defined $test{"$f[0] @f[4..6] $f[9] $f[16]"}); $test{"$f[0] @f[4..6] $f[9] $f[16]"} = 1;' > dia+rel/rem_dia_rel.maf

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