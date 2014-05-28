export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work

all: dia rel ../gene-patient-matrix.annotated.tsv ../gene-patient-matrix.af20.annotated.tsv pathscan rel/mutation-relation.tsv rel/mutation-relation.KRAS-PTPN11.tsv  

dia: ucsc-genes.hg19.roi dia/rem_dia.maf dia/gene_mrs dia/smg.tsv pathscan-dia
rel: ucsc-genes.hg19.roi rel/rem_rel.maf rel/gene_mrs rel/smg.tsv pathscan-rel
dia-af20: ucsc-genes.hg19.roi dia-af20/rem_dia.maf dia-af20/gene_mrs dia-af20/smg.tsv pathscan-dia
rel-af20: ucsc-genes.hg19.roi rel-af20/rem_rel.maf rel-af20/gene_mrs rel-af20/smg.tsv pathscan-rel
dia-af20-del: ucsc-genes.hg19.roi dia-af20-del/rem_dia.maf dia-af20-del/gene_mrs dia-af20-del/smg.tsv pathscan-dia-del
rel-af20-del: ucsc-genes.hg19.roi rel-af20-del/rem_rel.maf rel-af20-del/gene_mrs rel-af20-del/smg.tsv pathscan-rel-del

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

#-----------------------------------------------------------------------
# PREPARE MAF FILES FOR GENOME MUSIC
#-----------------------------------------------------------------------

PATIENTS = 314 1021247 399 430 446 460 545 592 715 786 792 818 842 A B C D E X Y

dia/rem_dia.maf: $(foreach P, $(PATIENTS), maf/$P_rem_dia.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>dia/rem_dia.maf.part | tee -a make.log
	cat maf/*rem_dia.maf >> dia/rem_dia.maf.part
	mv dia/rem_dia.maf.part dia/rem_dia.maf

rel/rem_rel.maf: $(foreach P, $(PATIENTS), maf/$P_rem_rel.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>rel/rem_rel.maf.part | tee -a make.log
	cat maf/*rem_rel.maf >> rel/rem_rel.maf.part
	mv rel/rem_rel.maf.part rel/rem_rel.maf

dia-af20/rem_dia.maf: $(foreach P, $(PATIENTS), maf-af20/$P_rem_dia.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>dia-af20/rem_dia.maf.part | tee -a make.log
	cat maf-af20/*rem_dia.maf >> dia-af20/rem_dia.maf.part
	mv dia-af20/rem_dia.maf.part dia-af20/rem_dia.maf

dia-af20/rem_dia.cnv.maf: dia-af20/rem_dia.maf ../cnv/hdall.cnv.tsv ~/hdall/scripts/pathway-analysis/cnv2maf.pl
	cp dia-af20/rem_dia.maf dia-af20/rem_dia.cnv.maf.part
	cat ../cnv/hdall.cnv.tsv | perl ~/hdall/scripts/pathway-analysis/cnv2maf.pl --max-size-bp 10000000 | grep "_dia" | grep "del" >> dia-af20/rem_dia.cnv.maf.part
	mv dia-af20/rem_dia.cnv.maf.part dia-af20/rem_dia.cnv.maf

dia-af20-del/rem_dia.maf: $(foreach P, $(PATIENTS), maf-af20-del/$P_rem_dia.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>dia-af20-del/rem_dia.maf.part | tee -a make.log
	cat maf-af20-del/*rem_dia.maf >> dia-af20-del/rem_dia.maf.part
	mv dia-af20-del/rem_dia.maf.part dia-af20-del/rem_dia.maf

dia-af20-del/rem_dia.cnv.maf: dia-af20-del/rem_dia.maf ../cnv/hdall.cnv.tsv ~/hdall/scripts/pathway-analysis/cnv2maf.pl
	cp dia-af20-del/rem_dia.maf dia-af20-del/rem_dia.cnv.maf.part
	cat ../cnv/hdall.cnv.tsv | perl ~/hdall/scripts/pathway-analysis/cnv2maf.pl --max-size-bp 10000000 | grep "_dia" | grep "del" >> dia-af20-del/rem_dia.cnv.maf.part
	mv dia-af20-del/rem_dia.cnv.maf.part dia-af20-del/rem_dia.cnv.maf

rel-af20/rem_rel.maf: $(foreach P, $(PATIENTS), maf-af20/$P_rem_rel.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>rel-af20/rem_rel.maf.part | tee -a make.log
	cat maf-af20/*rem_rel.maf >> rel-af20/rem_rel.maf.part
	mv rel-af20/rem_rel.maf.part rel-af20/rem_rel.maf

rel-af20/rem_rel.cnv.maf: rel-af20/rem_rel.maf ../cnv/hdall.cnv.tsv ~/hdall/scripts/pathway-analysis/cnv2maf.pl
	cp rel-af20/rem_rel.maf rel-af20/rem_rel.cnv.maf.part
	cat ../cnv/hdall.cnv.tsv | perl ~/hdall/scripts/pathway-analysis/cnv2maf.pl --max-size-bp 10000000 | grep "_rel" | grep "del" >> rel-af20/rem_rel.cnv.maf.part
	mv rel-af20/rem_rel.cnv.maf.part rel-af20/rem_rel.cnv.maf

rel-af20-del/rem_rel.maf: $(foreach P, $(PATIENTS), maf-af20-del/$P_rem_rel.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>rel-af20-del/rem_rel.maf.part | tee -a make.log
	cat maf-af20-del/*rem_rel.maf >> rel-af20-del/rem_rel.maf.part
	mv rel-af20-del/rem_rel.maf.part rel-af20-del/rem_rel.maf

rel-af20-del/rem_rel.cnv.maf: rel-af20-del/rem_rel.maf ../cnv/hdall.cnv.tsv ~/hdall/scripts/pathway-analysis/cnv2maf.pl
	cp rel-af20-del/rem_rel.maf rel-af20-del/rem_rel.cnv.maf.part
	cat ../cnv/hdall.cnv.tsv | perl ~/hdall/scripts/pathway-analysis/cnv2maf.pl --max-size-bp 10000000 | grep "_rel" | grep "del" >> rel-af20-del/rem_rel.cnv.maf.part
	mv rel-af20-del/rem_rel.cnv.maf.part rel-af20-del/rem_rel.cnv.maf

maf/%.maf: ../filtered_variants/%.snp.filtered.vcf ../filtered_variants/%.indel.filtered.vcf ucsc-genes.hg19.roi ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv ~/hdall/scripts/vcf2maf.pl
	mkdir -p maf
	cat ../filtered_variants/$*.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		2>&1 1>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)' | tee -a make.log
	cat ../filtered_variants/$*.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		2>&1 1>>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)' | tee -a make.log
	mv $@.part $@

maf-af20/%.maf: ../filtered_variants/%.snp.filtered.vcf ../filtered_variants/%.indel.filtered.vcf ucsc-genes.hg19.roi ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv ~/hdall/scripts/vcf2maf.pl
	mkdir -p maf-af20
	cat ../filtered_variants/$*.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.20 \
		2>&1 1>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)' | tee -a make.log
	cat ../filtered_variants/$*.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.20 \
		2>&1 1>>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)' | tee -a make.log
	mv $@.part $@

maf-af20-del/%.maf: ../filtered_variants/%.snp.filtered.vcf ../filtered_variants/%.indel.filtered.vcf ucsc-genes.hg19.roi ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv ~/hdall/scripts/vcf2maf.pl
	mkdir -p maf-af20-del
	cat ../filtered_variants/$*.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.20 \
		--deleterious \
		2>&1 1>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)' | tee -a make.log
	cat ../filtered_variants/$*.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.20 \
		--deleterious \
		2>&1 1>>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)' | tee -a make.log
	mv $@.part $@

#-----------------------------------------------------------------------
# GENOME MUSIC CALC-WIG-COVG
#-----------------------------------------------------------------------
	
dia/total_covgs: dia/wig-list.tsv ~/generic/data/hg19/ucsc.hg19.nochr.fasta
	rm -f dia/total_covgs
	rm -rf dia/roi_covgs
	rm -rf dia/gene_covgs
	genome music bmr calc-wig-covg \
		--wig-list dia/wig-list.tsv \
		--output-dir dia/ \
		--reference-sequence ~/generic/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

rel/total_covgs: rel/wig-list.tsv ~/generic/data/hg19/ucsc.hg19.nochr.fasta
	rm -f rel/total_covgs
	rm -rf rel/roi_covgs
	rm -rf rel/gene_covgs
	genome music bmr calc-wig-covg \
		--wig-list rel/wig-list.tsv \
		--output-dir rel/ \
		--reference-sequence ~/generic/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

dia-af20/total_covgs: dia/total_covgs
	rm -f dia-af20/total_covgs
	rm -rf dia-af20/roi_covgs
	rm -rf dia-af20/gene_covgs
	cp dia/total_covgs dia-af20/total_covgs
	cp -r dia/roi_covgs dia-af20/roi_covgs
	cp -r dia/gene_covgs dia-af20/gene_covgs

dia-af20-del/total_covgs: dia/total_covgs
	rm -f dia-af20-del/total_covgs
	rm -rf dia-af20-del/roi_covgs
	rm -rf dia-af20-del/gene_covgs
	cp dia/total_covgs dia-af20-del/total_covgs
	cp -r dia/roi_covgs dia-af20-del/roi_covgs
	cp -r dia/gene_covgs dia-af20-del/gene_covgs

rel-af20/total_covgs: rel/total_covgs
	rm -f rel-af20/total_covgs
	rm -rf rel-af20/roi_covgs
	rm -rf rel-af20/gene_covgs
	cp rel/total_covgs rel-af20/total_covgs
	cp -r rel/roi_covgs rel-af20/roi_covgs
	cp -r rel/gene_covgs rel-af20/gene_covgs

rel-af20-del/total_covgs: rel/total_covgs
	rm -f rel-af20-del/total_covgs
	rm -rf rel-af20-del/roi_covgs
	rm -rf rel-af20-del/gene_covgs
	cp rel/total_covgs rel-af20-del/total_covgs
	cp -r rel/roi_covgs rel-af20-del/roi_covgs
	cp -r rel/gene_covgs rel-af20-del/gene_covgs

#-----------------------------------------------------------------------
# GENOME MUSIC CALC-BMR
#-----------------------------------------------------------------------

dia/gene_mrs: dia/total_covgs dia/wig-list.tsv dia/rem_dia.maf ~/generic/data/hg19/ucsc.hg19.nochr.fasta ucsc-genes.hg19.roi
	rm -f dia/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list dia/wig-list.tsv \
		--maf-file dia/rem_dia.maf \
		--output-dir dia/ \
		--reference-sequence ~/generic/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

rel/gene_mrs: rel/total_covgs rel/wig-list.tsv rel/rem_rel.maf ~/generic/data/hg19/ucsc.hg19.nochr.fasta ucsc-genes.hg19.roi
	rm -f rel/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list rel/wig-list.tsv \
		--maf-file rel/rem_rel.maf \
		--output-dir rel/ \
		--reference-sequence ~/generic/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

dia-af20/gene_mrs: dia-af20/total_covgs dia/wig-list.tsv dia-af20/rem_dia.maf ~/generic/data/hg19/ucsc.hg19.nochr.fasta ucsc-genes.hg19.roi
	rm -f dia-af20/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list dia/wig-list.tsv \
		--maf-file dia-af20/rem_dia.maf \
		--output-dir dia-af20/ \
		--reference-sequence ~/generic/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

dia-af20-del/gene_mrs: dia-af20-del/total_covgs dia/wig-list.tsv dia-af20-del/rem_dia.maf ~/generic/data/hg19/ucsc.hg19.nochr.fasta ucsc-genes.hg19.roi
	rm -f dia-af20-del/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list dia/wig-list.tsv \
		--maf-file dia-af20-del/rem_dia.maf \
		--output-dir dia-af20-del/ \
		--reference-sequence ~/generic/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

rel-af20/gene_mrs: rel-af20/total_covgs rel/wig-list.tsv rel-af20/rem_rel.maf ~/generic/data/hg19/ucsc.hg19.nochr.fasta ucsc-genes.hg19.roi
	rm -f rel-af20/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list rel/wig-list.tsv \
		--maf-file rel-af20/rem_rel.maf \
		--output-dir rel-af20/ \
		--reference-sequence ~/generic/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

rel-af20-del/gene_mrs: rel-af20-del/total_covgs rel/wig-list.tsv rel-af20-del/rem_rel.maf ~/generic/data/hg19/ucsc.hg19.nochr.fasta ucsc-genes.hg19.roi
	rm -f rel-af20-del/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list rel/wig-list.tsv \
		--maf-file rel-af20-del/rem_rel.maf \
		--output-dir rel-af20-del/ \
		--reference-sequence ~/generic/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

#-----------------------------------------------------------------------
#--- GENOME MUSIC PROXIMITY
#-----------------------------------------------------------------------

dia/proximity.tsv: dia/rem_dia.maf
	genome music proximity \
		--maf-file dia/rem_dia.maf \
		--max-proximity 15 \
		--output-dir dia/ \
		2>&1 | tee -a make.log
	mv dia/proximity_report dia/proximity.tsv

rel/proximity.tsv: rel/rem_rel.maf
	genome music proximity \
		--maf-file rel/rem_rel.maf \
		--max-proximity 15 \
		--output-dir rel/ \
		2>&1 | tee -a make.log	
	mv rel/proximity_report rel/proximity.tsv
	
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

dia-af20/smg.tsv: dia-af20/gene_mrs
	genome music smg \
		--gene-mr-file dia-af20/gene_mrs \
		--output-file dia-af20/smg.tsv.part
	mv dia-af20/smg.tsv.part dia-af20/smg.tsv

dia-af20-del/smg.tsv: dia-af20-del/gene_mrs
	genome music smg \
		--gene-mr-file dia-af20-del/gene_mrs \
		--output-file dia-af20-del/smg.tsv.part
	mv dia-af20-del/smg.tsv.part dia-af20-del/smg.tsv

rel-af20/smg.tsv: rel-af20/gene_mrs
	genome music smg \
		--gene-mr-file rel-af20/gene_mrs \
		--output-file rel-af20/smg.tsv.part
	mv rel-af20/smg.tsv.part rel-af20/smg.tsv

rel-af20-del/smg.tsv: rel-af20-del/gene_mrs
	genome music smg \
		--gene-mr-file rel-af20-del/gene_mrs \
		--output-file rel-af20-del/smg.tsv.part
	mv rel-af20-del/smg.tsv.part rel-af20-del/smg.tsv

#-----------------------------------------------------------------------
#--- GENOME MUSIC PATH-SCAN (using only high allelic frequency variants)
#-----------------------------------------------------------------------

pathscan-dia: pathscan/sm_pathways.dia.all.annotated.tsv pathscan/sm_pathways.dia.noRAS.noCREBBP.annotated.tsv pathscan/sm_pathways.dia.all.cnv.annotated.tsv pathscan/sm_pathways.dia.all.matrix.clustered.tsv pathscan/sm_pathways.dia.noRAS.matrix.clustered.tsv pathscan/sm_pathways.dia.noCREBBP.matrix.clustered.tsv pathscan/sm_pathways.dia.noRAS.noCREBBP.matrix.clustered.tsv  
pathscan-rel: pathscan/sm_pathways.rel.all.annotated.tsv pathscan/sm_pathways.rel.noRAS.noCREBBP.annotated.tsv pathscan/sm_pathways.rel.all.cnv.annotated.tsv pathscan/sm_pathways.rel.all.matrix.clustered.tsv pathscan/sm_pathways.rel.noRAS.matrix.clustered.tsv pathscan/sm_pathways.rel.noCREBBP.matrix.clustered.tsv pathscan/sm_pathways.rel.noRAS.noCREBBP.matrix.clustered.tsv
pathscan-dia-del: pathscan/sm_pathways.dia.all.del.matrix.clustered.tsv pathscan/sm_pathways.dia.all.del.annotated.tsv pathscan/sm_pathways.dia.all.cnv.del.annotated.tsv 
pathscan-rel-del: pathscan/sm_pathways.rel.all.del.matrix.clustered.tsv pathscan/sm_pathways.rel.all.del.annotated.tsv pathscan/sm_pathways.rel.all.cnv.del.annotated.tsv 
pathscan: pathscan-dia pathscan-rel pathscan-dia-del pathscan-rel-del pathscan/sm_pathways.dia+rel.all.tsv pathscan/sm_pathways.dia+rel.all.del.tsv pathscan/sm_pathways.dia+rel.noRAS.noCREBBP.tsv ../pathway-patient-matrix.all.clustered.tsv

pathscan/pathways-david-complete.tsv: pathscan/pathways-david.tsv pathscan/pathways-david-reverse.tsv
	perl ~/hdall/scripts/pathway-analysis/merge-david-annotations.pl \
		pathscan/pathways-david.tsv \
		pathscan/pathways-david-reverse.tsv \
		> pathscan/pathways-david-complete.tsv.part
	mv pathscan/pathways-david-complete.tsv.part pathscan/pathways-david-complete.tsv

pathscan/gene-sets.mutsigdb4.0.tsv: ~/generic/data/msigdb/msigdb.v4.0.entrez.gmt ~/hdall/scripts/pathway-analysis/msigdb2pathscan.pl
	cat ~/generic/data/msigdb/c1.all.v4.0.symbols.gmt | perl ~/hdall/scripts/pathway-analysis/msigdb2pathscan.pl \
		--entrez-id-file ~/generic/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c1_positional \
		2>&1 1>$@.part | tee -a make.log 
	cat ~/generic/data/msigdb/c2.all.v4.0.symbols.gmt | perl ~/hdall/scripts/pathway-analysis/msigdb2pathscan.pl \
		--entrez-id-file ~/generic/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c2_curated \
		2>&1 1>>$@.part | tee -a make.log 
	cat ~/generic/data/msigdb/c3.all.v4.0.symbols.gmt | perl ~/hdall/scripts/pathway-analysis/msigdb2pathscan.pl \
		--entrez-id-file ~/generic/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c3_motif \
		2>&1 1>>$@.part | tee -a make.log 
	cat ~/generic/data/msigdb/c4.all.v4.0.symbols.gmt | perl ~/hdall/scripts/pathway-analysis/msigdb2pathscan.pl \
		--entrez-id-file ~/generic/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c4_computational \
		2>&1 1>>$@.part | tee -a make.log 
	cat ~/generic/data/msigdb/c5.all.v4.0.symbols.gmt | perl ~/hdall/scripts/pathway-analysis/msigdb2pathscan.pl \
		--entrez-id-file ~/generic/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c5_GO \
		2>&1 1>>$@.part | tee -a make.log 
	cat ~/generic/data/msigdb/c6.all.v4.0.symbols.gmt | perl ~/hdall/scripts/pathway-analysis/msigdb2pathscan.pl \
		--entrez-id-file ~/generic/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c6_oncogenic_signature \
		2>&1 1>>$@.part | tee -a make.log 
	cat ~/generic/data/msigdb/c7.all.v4.0.symbols.gmt | perl ~/hdall/scripts/pathway-analysis/msigdb2pathscan.pl \
		--entrez-id-file ~/generic/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c7_immunologic_signature \
		2>&1 1>>$@.part | tee -a make.log
	mv $@.part $@ 

pathscan/pathways-nci.tsv: ~/generic/data/nci-pathway/uniprot.tab ~/hdall/scripts/pathway-analysis/nci2pathscan.pl
	cat ~/generic/data/nci-pathway/uniprot.tab | perl ~/hdall/scripts/pathway-analysis/nci2pathscan.pl \
		2>&1 1>>pathscan/pathways-nci.tsv.part | tee -a make.log
	mv pathscan/pathways-nci.tsv.part pathscan/pathways-nci.tsv

pathscan/pathways-wikipathways.tsv: ~/generic/data/wikipathways/wikipathways_Homo_sapiens_Curation-AnalysisCollection__txt.zip
	perl ~/hdall/scripts/pathway-analysis/wikipathways2pathscan.pl \
		--wiki-dir ~/generic/data/wikipathways/ \
		2>&1 1>>pathscan/pathways-wikipathways.tsv.part | tee -a make.log
	mv pathscan/pathways-wikipathways.tsv.part pathscan/pathways-wikipathways.tsv

pathscan/pathways-kegg.jerry.hugo.tsv: pathscan/pathways-kegg.jerry.ensembl.tsv
	perl ~/hdall/scripts/pathway-analysis/kegg-jerry-ensembl-to-hugo.pl > $@.part
	mv $@.part $@

pathscan/pathways-all.tsv: pathscan/pathways-david-complete.tsv pathscan/gene-sets.mutsigdb4.0.tsv pathscan/pathways-nci.tsv pathscan/pathways-wikipathways.tsv pathscan/pathways-kegg.jerry.hugo.tsv
	cp pathscan/pathways-david-complete.tsv pathscan/pathways-all.tsv.part
	cat pathscan/gene-sets.mutsigdb4.0.tsv pathscan/pathways-nci.tsv pathscan/pathways-manual.tsv pathscan/pathways-wikipathways.tsv >> pathscan/pathways-all.tsv.part
	cat pathscan/pathways-kegg.jerry.hugo.tsv >> pathscan/pathways-all.tsv.part
	mv pathscan/pathways-all.tsv.part pathscan/pathways-all.tsv

pathscan/pathways-all.sizefiltered.tsv: pathscan/pathways-all.tsv
	cat pathscan/pathways-all.tsv | perl -ne '@f = split /\t/; @g = split(/\|/, $$f[3]); print $$_ if (@g > 1 and @g < 400);' > pathscan/pathways-all.sizefiltered.tsv.part
	mv pathscan/pathways-all.sizefiltered.tsv.part pathscan/pathways-all.sizefiltered.tsv

pathscan/pathways-kegg-biocarta-bbid-reactome.tsv: pathscan/pathways-david-complete.tsv pathscan/gene-sets.mutsigdb4.0.tsv
	grep -P "(KEGG_PATHWAY|BIOCARTA|BBID)" pathscan/pathways-david-complete.tsv > pathscan/pathways-kegg-biocarta-bbid-reactome.tsv
	grep "REACTOME_" pathscan/gene-sets.mutsigdb4.0.tsv >> pathscan/pathways-kegg-biocarta-bbid-reactome.tsv

# NOTE: only exclude likely false positives from Lawrence et al. (2013)
pathscan/sm_pathways.dia.all.raw.tsv: dia-af20/smg.tsv pathscan/pathways-all.sizefiltered.tsv dia-af20/rem_dia.maf
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia-af20/gene_covgs/ \
		--maf-file dia-af20/rem_dia.maf \
		--output-file pathscan/sm_pathways.dia.all.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.dia.all.raw.tsv.part pathscan/sm_pathways.dia.all.raw.tsv
	mv pathscan/sm_pathways.dia.all.raw.tsv.part_detailed pathscan/sm_pathways.dia.all.raw.tsv_detailed

pathscan/sm_pathways.dia.all.del.raw.tsv: dia-af20-del/smg.tsv pathscan/pathways-all.sizefiltered.tsv dia-af20-del/rem_dia.maf
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia-af20-del/gene_covgs/ \
		--maf-file dia-af20-del/rem_dia.maf \
		--output-file pathscan/sm_pathways.dia.all.del.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.dia.all.del.raw.tsv.part pathscan/sm_pathways.dia.all.del.raw.tsv
	mv pathscan/sm_pathways.dia.all.del.raw.tsv.part_detailed pathscan/sm_pathways.dia.all.del.raw.tsv_detailed

pathscan/sm_pathways.rel.all.raw.tsv: rel-af20/smg.tsv pathscan/pathways-all.sizefiltered.tsv rel-af20/rem_rel.maf
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel-af20/gene_covgs/ \
		--maf-file rel-af20/rem_rel.maf \
		--output-file pathscan/sm_pathways.rel.all.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.rel.all.raw.tsv.part pathscan/sm_pathways.rel.all.raw.tsv
	mv pathscan/sm_pathways.rel.all.raw.tsv.part_detailed pathscan/sm_pathways.rel.all.raw.tsv_detailed

pathscan/sm_pathways.rel.all.del.raw.tsv: rel-af20-del/smg.tsv pathscan/pathways-all.sizefiltered.tsv rel-af20-del/rem_rel.maf
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel-af20-del/gene_covgs/ \
		--maf-file rel-af20-del/rem_rel.maf \
		--output-file pathscan/sm_pathways.rel.all.del.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.rel.all.del.raw.tsv.part pathscan/sm_pathways.rel.all.del.raw.tsv
	mv pathscan/sm_pathways.rel.all.del.raw.tsv.part_detailed pathscan/sm_pathways.rel.all.del.raw.tsv_detailed

# NOTE: consider also CNV data from arrays
pathscan/sm_pathways.dia.all.cnv.raw.tsv: dia-af20/smg.tsv pathscan/pathways-all.sizefiltered.tsv dia-af20/rem_dia.cnv.maf
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia-af20/gene_covgs/ \
		--maf-file dia-af20/rem_dia.cnv.maf \
		--output-file pathscan/sm_pathways.dia.all.cnv.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.dia.all.cnv.raw.tsv.part pathscan/sm_pathways.dia.all.cnv.raw.tsv
	mv pathscan/sm_pathways.dia.all.cnv.raw.tsv.part_detailed pathscan/sm_pathways.dia.all.cnv.raw.tsv_detailed

pathscan/sm_pathways.dia.all.cnv.del.raw.tsv: dia-af20-del/smg.tsv pathscan/pathways-all.sizefiltered.tsv dia-af20-del/rem_dia.cnv.maf
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia-af20-del/gene_covgs/ \
		--maf-file dia-af20-del/rem_dia.cnv.maf \
		--output-file pathscan/sm_pathways.dia.all.cnv.del.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.dia.all.cnv.del.raw.tsv.part pathscan/sm_pathways.dia.all.cnv.del.raw.tsv
	mv pathscan/sm_pathways.dia.all.cnv.del.raw.tsv.part_detailed pathscan/sm_pathways.dia.all.cnv.del.raw.tsv_detailed

pathscan/sm_pathways.rel.all.cnv.raw.tsv: rel-af20/smg.tsv pathscan/pathways-all.sizefiltered.tsv rel-af20/rem_rel.cnv.maf
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel-af20/gene_covgs/ \
		--maf-file rel-af20/rem_rel.cnv.maf \
		--output-file pathscan/sm_pathways.rel.all.cnv.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.rel.all.cnv.raw.tsv.part pathscan/sm_pathways.rel.all.cnv.raw.tsv
	mv pathscan/sm_pathways.rel.all.cnv.raw.tsv.part_detailed pathscan/sm_pathways.rel.all.cnv.raw.tsv_detailed

pathscan/sm_pathways.rel.all.cnv.del.raw.tsv: rel-af20-del/smg.tsv pathscan/pathways-all.sizefiltered.tsv rel-af20-del/rem_rel.cnv.maf
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel-af20-del/gene_covgs/ \
		--maf-file rel-af20-del/rem_rel.cnv.maf \
		--output-file pathscan/sm_pathways.rel.all.cnv.del.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.rel.all.cnv.del.raw.tsv.part pathscan/sm_pathways.rel.all.cnv.del.raw.tsv
	mv pathscan/sm_pathways.rel.all.cnv.del.raw.tsv.part_detailed pathscan/sm_pathways.rel.all.cnv.del.raw.tsv_detailed

# NOTE: exclude KRAS,NRAS
pathscan/sm_pathways.dia.noRAS.raw.tsv: dia-af20-del/smg.tsv pathscan/pathways-all.sizefiltered.tsv dia-af20-del/rem_dia.maf
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia-af20-del/gene_covgs/ \
		--maf-file dia-af20-del/rem_dia.maf \
		--output-file pathscan/sm_pathways.dia.noRAS.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,KRAS,NRAS,PTPN11,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.dia.noRAS.raw.tsv.part pathscan/sm_pathways.dia.noRAS.raw.tsv
	mv pathscan/sm_pathways.dia.noRAS.raw.tsv.part_detailed pathscan/sm_pathways.dia.noRAS.raw.tsv_detailed

pathscan/sm_pathways.rel.noRAS.raw.tsv: rel-af20-del/smg.tsv pathscan/pathways-all.sizefiltered.tsv rel-af20-del/rem_rel.maf
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel-af20-del/gene_covgs/ \
		--maf-file rel-af20-del/rem_rel.maf \
		--output-file pathscan/sm_pathways.rel.noRAS.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,KRAS,NRAS,PTPN11,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.rel.noRAS.raw.tsv.part pathscan/sm_pathways.rel.noRAS.raw.tsv
	mv pathscan/sm_pathways.rel.noRAS.raw.tsv.part_detailed pathscan/sm_pathways.rel.noRAS.raw.tsv_detailed

# NOTE: exclude CREBBP
pathscan/sm_pathways.dia.noCREBBP.raw.tsv: dia-af20-del/smg.tsv pathscan/pathways-all.sizefiltered.tsv dia-af20-del/rem_dia.maf
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia-af20-del/gene_covgs/ \
		--maf-file dia-af20-del/rem_dia.maf \
		--output-file pathscan/sm_pathways.dia.noCREBBP.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,CREBBP,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.dia.noCREBBP.raw.tsv.part pathscan/sm_pathways.dia.noCREBBP.raw.tsv
	mv pathscan/sm_pathways.dia.noCREBBP.raw.tsv.part_detailed pathscan/sm_pathways.dia.noCREBBP.raw.tsv_detailed

pathscan/sm_pathways.rel.noCREBBP.raw.tsv: rel-af20-del/smg.tsv pathscan/pathways-all.sizefiltered.tsv rel-af20-del/rem_rel.maf
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel-af20-del/gene_covgs/ \
		--maf-file rel-af20-del/rem_rel.maf \
		--output-file pathscan/sm_pathways.rel.noCREBBP.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,CREBBP,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.rel.noCREBBP.raw.tsv.part pathscan/sm_pathways.rel.noCREBBP.raw.tsv
	mv pathscan/sm_pathways.rel.noCREBBP.raw.tsv.part_detailed pathscan/sm_pathways.rel.noCREBBP.raw.tsv_detailed

# NOTE: exclude KRAS,NRAS,CREBBP
pathscan/sm_pathways.dia.noRAS.noCREBBP.raw.tsv: dia-af20-del/smg.tsv pathscan/pathways-all.sizefiltered.tsv dia-af20-del/rem_dia.maf
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia-af20-del/gene_covgs/ \
		--maf-file dia-af20-del/rem_dia.maf \
		--output-file pathscan/sm_pathways.dia.noRAS.noCREBBP.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,KRAS,NRAS,PTPN11,CREBBP,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.dia.noRAS.noCREBBP.raw.tsv.part pathscan/sm_pathways.dia.noRAS.noCREBBP.raw.tsv
	mv pathscan/sm_pathways.dia.noRAS.noCREBBP.raw.tsv.part_detailed pathscan/sm_pathways.dia.noRAS.noCREBBP.raw.tsv_detailed

pathscan/sm_pathways.rel.noRAS.noCREBBP.raw.tsv: rel-af20-del/smg.tsv pathscan/pathways-all.sizefiltered.tsv rel-af20-del/rem_rel.maf
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel-af20-del/gene_covgs/ \
		--maf-file rel-af20-del/rem_rel.maf \
		--output-file pathscan/sm_pathways.rel.noRAS.noCREBBP.raw.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore MESP2,KRAS,NRAS,PTPN11,CREBBP,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1,PDZD7
	mv pathscan/sm_pathways.rel.noRAS.noCREBBP.raw.tsv.part pathscan/sm_pathways.rel.noRAS.noCREBBP.raw.tsv
	mv pathscan/sm_pathways.rel.noRAS.noCREBBP.raw.tsv.part_detailed pathscan/sm_pathways.rel.noRAS.noCREBBP.raw.tsv_detailed

pathscan/sm_pathways.dia.%.annotated.tsv: pathscan/sm_pathways.dia.%.raw.tsv ~/hdall/scripts/pathway-analysis/annotate-pathscan-result.pl
	perl ~/hdall/scripts/pathway-analysis/annotate-pathscan-result.pl \
		--pathway-file pathscan/pathways-all.tsv \
		--sm-pathways pathscan/sm_pathways.dia.$*.raw.tsv \
		--sm-pathways-detail pathscan/sm_pathways.dia.$*.raw.tsv_detailed \
		--maf-file dia-af20/rem_dia.cnv.maf \
		> pathscan/sm_pathways.dia.$*.annotated.tsv.part
	mv pathscan/sm_pathways.dia.$*.annotated.tsv.part pathscan/sm_pathways.dia.$*.annotated.tsv

pathscan/sm_pathways.rel.%.annotated.tsv: pathscan/sm_pathways.rel.%.raw.tsv ~/hdall/scripts/pathway-analysis/annotate-pathscan-result.pl
	perl ~/hdall/scripts/pathway-analysis/annotate-pathscan-result.pl \
		--pathway-file pathscan/pathways-all.tsv \
		--sm-pathways pathscan/sm_pathways.rel.$*.raw.tsv \
		--sm-pathways-detail pathscan/sm_pathways.rel.$*.raw.tsv_detailed \
		--maf-file rel-af20/rem_rel.cnv.maf \
		> pathscan/sm_pathways.rel.$*.annotated.tsv.part
	mv pathscan/sm_pathways.rel.$*.annotated.tsv.part pathscan/sm_pathways.rel.$*.annotated.tsv

pathscan/sm_pathways.%.matrix.clustered.tsv: pathscan/sm_pathways.%.annotated.tsv ~/hdall/scripts/pathway-analysis/get-gene-pathway-matrix-enriched.pl ~/hdall/scripts/pathway-analysis/cluster-gene-pathway-matrix.R
	cat pathscan/sm_pathways.$*.annotated.tsv | perl ~/hdall/scripts/pathway-analysis/get-gene-pathway-matrix-enriched.pl \
		--max-p-value 1e-6 \
		--max-pathway-size 400 \
		2>&1 1>pathscan/sm_pathways.$*.matrix.tsv.tmp | tee -a make.log
	rm -f pathscan/sm_pathways.$*.matrix.genetree.pdf pathscan/sm_pathways.$*.matrix.pathwaytree.pdf
	R --no-save --quiet --slave -f ~/hdall/scripts/pathway-analysis/cluster-gene-pathway-matrix.R \
		--args pathscan/sm_pathways.$*.matrix.tsv.tmp pathscan/sm_pathways.$*.matrix.clustered.tsv.part pathscan/sm_pathways.$*.matrix.genetree.pdf pathscan/sm_pathways.$*.matrix.pathwaytree.pdf \
		2>&1 | tee -a make.log
	rm pathscan/sm_pathways.$*.matrix.tsv.tmp
	mv pathscan/sm_pathways.$*.matrix.clustered.tsv.part pathscan/sm_pathways.$*.matrix.clustered.tsv

pathscan/sm_pathways.rel-not-dia.all.matrix.clustered.tsv: pathscan/sm_pathways.dia.all.annotated.tsv pathscan/sm_pathways.rel.all.annotated.tsv ~/hdall/scripts/pathway-analysis/get-gene-pathway-matrix-enriched.pl ~/hdall/scripts/pathway-analysis/cluster-gene-pathway-matrix.R
	cat pathscan/sm_pathways.rel.all.annotated.tsv | perl ~/hdall/scripts/pathway-analysis/get-gene-pathway-matrix-enriched.pl \
		--max-p-value 1e-5 \
		--max-pathway-size 400 \
		--subtract-pathways pathscan/sm_pathways.dia.all.annotated.tsv \
		2>&1 1>pathscan/sm_pathways.rel-not-dia.all.matrix.tsv.tmp | tee -a make.log
	R --no-save --quiet --slave -f ~/hdall/scripts/pathway-analysis/cluster-gene-pathway-matrix.R \
		--args pathscan/sm_pathways.rel-not-dia.all.matrix.tsv.tmp pathscan/sm_pathways.rel-not-dia.all.matrix.clustered.tsv pathscan/sm_pathways.rel-not-dia.all.matrix.genetree.pdf pathscan/sm_pathways.rel-not-dia.all.matrix.pathwaytree.pdf \
		2>&1 | tee -a make.log
	rm pathscan/sm_pathways.rel-not-dia.all.matrix.tsv.tmp

pathscan/sm_pathways.dia+rel.all.tsv: pathscan/sm_pathways.dia.all.annotated.tsv pathscan/sm_pathways.rel.all.annotated.tsv ~/hdall/scripts/pathway-analysis/cluster-pathways.R ~/hdall/scripts/music/merge-sm-pathways.R
	R --no-save --quiet --slave -f ~/hdall/scripts/music/merge-sm-pathways.R \
		--args pathscan/sm_pathways.dia.all.annotated.tsv pathscan/sm_pathways.rel.all.annotated.tsv \
		2>&1 1>pathscan/sm_pathways.dia+rel.all.tsv.part | tee -a make.log
	R --no-save --quiet --slave -f ~/hdall/scripts/pathway-analysis/cluster-pathways.R \
		--args pathscan/sm_pathways.dia+rel.all.tsv.part pathscan/sm_pathways.dia+rel.all.clustered.tsv.part 1e-7 \
		2>&1 | tee -a make.log
	mv pathscan/sm_pathways.dia+rel.all.clustered.tsv.part pathscan/sm_pathways.dia+rel.all.tsv
	rm pathscan/sm_pathways.dia+rel.all.tsv.part

pathscan/sm_pathways.dia+rel.all.del.tsv: pathscan/sm_pathways.dia.all.del.annotated.tsv pathscan/sm_pathways.rel.all.del.annotated.tsv ~/hdall/scripts/pathway-analysis/cluster-pathways.R ~/hdall/scripts/music/merge-sm-pathways.R
	R --no-save --quiet --slave -f ~/hdall/scripts/music/merge-sm-pathways.R \
		--args pathscan/sm_pathways.dia.all.del.annotated.tsv pathscan/sm_pathways.rel.all.del.annotated.tsv \
		2>&1 1>$@.part | tee -a make.log
	R --no-save --quiet --slave -f ~/hdall/scripts/pathway-analysis/cluster-pathways.R \
		--args $@.part pathscan/sm_pathways.dia+rel.all.del.clustered.tsv.part 1e-7 \
		2>&1 | tee -a make.log
	mv pathscan/sm_pathways.dia+rel.all.del.clustered.tsv.part $@
	rm $@.part

pathscan/sm_pathways.dia+rel.noRAS.noCREBBP.tsv: pathscan/sm_pathways.dia.noRAS.noCREBBP.annotated.tsv pathscan/sm_pathways.rel.noRAS.noCREBBP.annotated.tsv ~/hdall/scripts/pathway-analysis/cluster-pathways.R ~/hdall/scripts/music/merge-sm-pathways.R
	R --no-save --quiet --slave -f ~/hdall/scripts/music/merge-sm-pathways.R \
		--args pathscan/sm_pathways.dia.noRAS.noCREBBP.annotated.tsv pathscan/sm_pathways.rel.noRAS.noCREBBP.annotated.tsv \
		2>&1 1>$@.part | tee -a make.log
	R --no-save --quiet --slave -f ~/hdall/scripts/pathway-analysis/cluster-pathways.R \
		--args $@.part pathscan/sm_pathways.dia+rel.noRAS.noCREBBP.clustered.tsv.part 1e-7 \
		2>&1 | tee -a make.log
	mv pathscan/sm_pathways.dia+rel.noRAS.noCREBBP.clustered.tsv.part $@
	rm $@.part

../pathway-patient-matrix.all.clustered.tsv: pathscan/sm_pathways.dia.all.annotated.tsv pathscan/sm_pathways.rel.all.annotated.tsv ../filtered-variants.tsv ../kamilla/gene-pathway-assignment.tsv ~/hdall/scripts/pathway-analysis/get-pathway-patient-matrix.pl ~/hdall/scripts/pathway-analysis/cluster-pathways.R
	perl ~/hdall/scripts/pathway-analysis/get-pathway-patient-matrix.pl \
		--enriched-pathways-dia pathscan/sm_pathways.dia.all.annotated.tsv \
		--enriched-pathways-rel pathscan/sm_pathways.rel.all.annotated.tsv \
		--curated-pathways ../kamilla/gene-pathway-assignment.tsv \
		--filtered-variants ../filtered-variants.tsv \
		2>&1 1>../pathway-patient-matrix.tsv.part | tee -a make.log
	mv ../pathway-patient-matrix.tsv.part ../pathway-patient-matrix.tsv
	# TABLE: pathway-patient-matrix
	R --no-save --quiet --slave -f ~/hdall/scripts/pathway-analysis/cluster-pathways.R \
		--args ../pathway-patient-matrix.tsv ../pathway-patient-matrix.all.clustered.tsv.part 1e-7 \
		2>&1 | tee -a make.log
	mv ../pathway-patient-matrix.all.clustered.tsv.part ../pathway-patient-matrix.all.clustered.tsv

#-----------------------------------------------------------------------
#--- MUSIC MUTATION-RELATION ANALYSIS
#-----------------------------------------------------------------------

rel/mutation-relation.tsv: ../gene-patient-matrix.tsv rel/wig-list.tsv rel/rem_rel.maf
	rm -f rel/mutation-relation.gene-list-recrel2-maxaf20.txt
	python ~/hdall/scripts/query-text.py \
		--file ../gene-patient-matrix.tsv \
		--sql "select gene,data.'max-af-rel-ns' AS m,data.'freq-rel-ns' AS f from data where m >= 20 and m <> '' and f >= 2" \
		> rel/mutation-relation.gene-list-recrel2-maxaf20.txt
	genome music mutation-relation \
		--bam-list=rel/wig-list.tsv \
		--maf-file=rel/rem_rel.maf \
		--permutations 1000 \
		--gene-list rel/mutation-relation.gene-list-recrel2-maxaf20.txt \
		--output-file rel/mutation-relation.tsv.part
	mv rel/mutation-relation.tsv.part rel/mutation-relation.tsv

rel/mutation-relation.KRAS-PTPN11.tsv: ../gene-patient-matrix.tsv rel/wig-list.tsv rel/rem_rel.KRAS-PTPN11.maf rel/mutation-relation.gene-list-recrel2-maxaf20.KRAS-PTPN11.txt
	genome music mutation-relation \
		--bam-list=rel/wig-list.tsv \
		--maf-file=rel/rem_rel.KRAS-PTPN11.maf \
		--permutations 1000 \
		--gene-list rel/mutation-relation.gene-list-recrel2-maxaf20.KRAS-PTPN11.txt \
		--output-file rel/mutation-relation.KRAS-PTPN11.tsv.part
	mv rel/mutation-relation.KRAS-PTPN11.tsv.part rel/mutation-relation.KRAS-PTPN11.tsv

dia/mutation-relation.tsv: ../gene-patient-matrix.tsv dia/wig-list.tsv dia/rem_dia.maf
	rm -f dia/mutation-relation.gene-list-recdia2-maxaf20.txt
	python ~/hdall/scripts/query-text.py \
		--file ../gene-patient-matrix.tsv \
		--sql "select gene,data.'max-af-dia-ns' AS m,data.'freq-dia-ns' AS f from data where m >= 20 and m <> '' and f >= 2" \
		> dia/mutation-relation.gene-list-recdia2-maxaf20.txt
	genome music mutation-relation \
		--bam-list=dia/wig-list.tsv \
		--maf-file=dia/rem_dia.maf \
		--permutations 1000 \
		--gene-list dia/mutation-relation.gene-list-recdia2-maxaf20.txt \
		--output-file dia/mutation-relation.tsv.part
	mv dia/mutation-relation.tsv.part dia/mutation-relation.tsv

#-----------------------------------------------------------------------
#--- COMBINE ALL RESULTS INTO FINAL LISTS
#-----------------------------------------------------------------------

#pathscan/gene-pathway-matrix.clustered.tsv: pathscan/pathways-nci.tsv pathscan/pathways-kegg-biocarta-bbid-reactome.tsv pathscan/pathways-wikipathways.tsv ~/hdall/scripts/pathway-analysis/get-gene-pathway-matrix-full.pl ~/hdall/scripts/pathway-analysis/cluster-gene-pathway-matrix.R 
#	cat pathscan/pathways-nci.tsv pathscan/pathways-kegg-biocarta-bbid-reactome.tsv pathscan/pathways-wikipathways.tsv \
#		| perl ~/hdall/scripts/pathway-analysis/get-gene-pathway-matrix-full.pl --pathway-file - > pathscan/gene-pathway-matrix.tsv.part
#	mv pathscan/gene-pathway-matrix.tsv.part pathscan/gene-pathway-matrix.tsv
#	R --no-save --quiet --slave -f ~/hdall/scripts/pathway-analysis/cluster-gene-pathway-matrix.R \
#		--args pathscan/gene-pathway-matrix.tsv pathscan/gene-pathway-matrix.clustered.tsv \
#		2>&1 | tee -a make.log

../gene-patient-matrix.annotated.tsv: ../gene-patient-matrix.tsv dia-af20/smg.tsv rel-af20/smg.tsv pathscan/sm_pathways.dia.all.annotated.tsv pathscan/sm_pathways.rel.all.annotated.tsv ../cnv/impacted-genes-list.cnv.tsv pathscan/sm_pathways.rel.all.matrix.clustered.tsv ../freq-rel-ns-curated.tsv ~/hdall/scripts/annotate-gene-patient-matrix.pl
	perl ~/hdall/scripts/annotate-gene-patient-matrix.pl \
		--gene-patient-matrix ../gene-patient-matrix.tsv \
		--smg-dia dia-af20/smg.tsv \
		--smg-rel rel-af20/smg.tsv \
		--smp-dia pathscan/sm_pathways.dia.all.annotated.tsv \
		--smp-rel pathscan/sm_pathways.rel.all.annotated.tsv \
		--cnv ../cnv/impacted-genes-list.cnv.tsv \
		--gene-pathway-matrix pathscan/sm_pathways.rel.all.matrix.clustered.tsv \
		--panel-genes ../panel-genes.tsv \
		--max-p-value 1e-8 \
		--max-pathway-size 400 \
		--curated-frequencies ../freq-rel-ns-curated.tsv \
		2>&1 1>../gene-patient-matrix.annotated.tsv.part | tee -a make.log
	mv ../gene-patient-matrix.annotated.tsv.part ../gene-patient-matrix.annotated.tsv	
	
../gene-patient-matrix.af20.annotated.tsv: ../gene-patient-matrix.af20.tsv dia-af20/smg.tsv rel-af20/smg.tsv pathscan/sm_pathways.dia.noRAS.noCREBBP.annotated.tsv pathscan/sm_pathways.rel.noRAS.noCREBBP.annotated.tsv ../cnv/impacted-genes-list.cnv.tsv pathscan/sm_pathways.rel.noRAS.noCREBBP.matrix.clustered.tsv ~/hdall/scripts/annotate-gene-patient-matrix.pl
	perl ~/hdall/scripts/annotate-gene-patient-matrix.pl \
		--gene-patient-matrix ../gene-patient-matrix.af20.tsv \
		--smg-dia dia-af20-del/smg.tsv \
		--smg-rel rel-af20-del/smg.tsv \
		--smp-dia pathscan/sm_pathways.dia.noRAS.noCREBBP.annotated.tsv \
		--smp-rel pathscan/sm_pathways.rel.noRAS.noCREBBP.annotated.tsv \
		--cnv ../cnv/impacted-genes-list.cnv.tsv \
		--gene-pathway-matrix pathscan/sm_pathways.rel.noRAS.noCREBBP.matrix.clustered.tsv \
		--panel-genes ../panel-genes.tsv \
		--max-p-value 1e-5 \
		--max-pathway-size 400 \
		2>&1 1>../gene-patient-matrix.af20.annotated.tsv.part | tee -a make.log
	mv ../gene-patient-matrix.af20.annotated.tsv.part ../gene-patient-matrix.af20.annotated.tsv
