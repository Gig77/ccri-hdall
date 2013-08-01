export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work

all: dia rel ../gene-patient-matrix.annotated.tsv ../gene-patient-matrix.high-af.annotated.tsv pathscan/sm_pathways.tsv rel/mutation-relation.tsv rel/mutation-relation.KRAS-PTPN11.tsv ../pathway-patient-matrix.clustered.tsv ../pathway-patient-matrix.iter2.clustered.tsv  

dia: ucsc-genes.hg19.roi dia/rem_dia.maf dia/gene_mrs dia/smg.tsv pathscan/sm_pathways.dia.tsv pathscan/sm_pathways.dia.annotated.tsv
rel: ucsc-genes.hg19.roi rel/rem_rel.maf rel/gene_mrs rel/smg.tsv pathscan/sm_pathways.rel.tsv pathscan/sm_pathways.rel.annotated.tsv
dia-high-af: ucsc-genes.hg19.roi dia-high-af/rem_dia.maf dia-high-af/gene_mrs dia-high-af/smg.tsv pathscan/sm_pathways.dia.tsv pathscan/sm_pathways.dia.annotated.tsv
rel-high-af: ucsc-genes.hg19.roi rel-high-af/rem_rel.maf rel-high-af/gene_mrs rel-high-af/smg.tsv pathscan/sm_pathways.rel.tsv pathscan/sm_pathways.rel.annotated.tsv

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

dia-high-af/rem_dia.maf: $(foreach P, $(PATIENTS), maf-high-af/$P_rem_dia.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>dia-high-af/rem_dia.maf.part | tee -a make.log
	cat maf-high-af/*rem_dia.maf >> dia-high-af/rem_dia.maf.part
	mv dia-high-af/rem_dia.maf.part dia-high-af/rem_dia.maf

rel-high-af/rem_rel.maf: $(foreach P, $(PATIENTS), maf-high-af/$P_rem_rel.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>rel-high-af/rem_rel.maf.part | tee -a make.log
	cat maf-high-af/*rem_rel.maf >> rel-high-af/rem_rel.maf.part
	mv rel-high-af/rem_rel.maf.part rel-high-af/rem_rel.maf

maf/%.maf: ../filtered_variants/%.snp.filtered.vcf ../filtered_variants/%.indel.filtered.vcf ucsc-genes.hg19.roi ~/hdall/results/gene-id-mapping.biomart-0.7.tsv ~/hdall/scripts/vcf2maf.pl
	mkdir -p maf
	cat ../filtered_variants/$*.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		2>&1 1>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)' | tee -a make.log
	cat ../filtered_variants/$*.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		2>&1 1>>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)' | tee -a make.log
	mv $@.part $@

maf-high-af/%.maf: ../filtered_variants/%.snp.filtered.vcf ../filtered_variants/%.indel.filtered.vcf ucsc-genes.hg19.roi ~/hdall/results/gene-id-mapping.biomart-0.7.tsv ~/hdall/scripts/vcf2maf.pl
	mkdir -p maf-high-af
	cat ../filtered_variants/$*.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.20 \
		2>&1 1>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)' | tee -a make.log
	cat ../filtered_variants/$*.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.20 \
		2>&1 1>>$@.part | grep -vP '(Leading or trailing space|variant.Format|Domain annotations come from)' | tee -a make.log
	mv $@.part $@
		
#-----------------------------------------------------------------------
# GENOME MUSIC CALC-WIG-COVG
#-----------------------------------------------------------------------
	
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

dia-high-af/total_covgs: dia/total_covgs
	rm -f dia-high-af/total_covgs
	rm -rf dia-high-af/roi_covgs
	rm -rf dia-high-af/gene_covgs
	cp dia/total_covgs dia-high-af/total_covgs
	cp -r dia/roi_covgs dia-high-af/roi_covgs
	cp -r dia/gene_covgs dia-high-af/gene_covgs

rel-high-af/total_covgs: rel/total_covgs
	rm -f rel-high-af/total_covgs
	rm -rf rel-high-af/roi_covgs
	rm -rf rel-high-af/gene_covgs
	cp rel/total_covgs rel-high-af/total_covgs
	cp -r rel/roi_covgs rel-high-af/roi_covgs
	cp -r rel/gene_covgs rel-high-af/gene_covgs

#-----------------------------------------------------------------------
# GENOME MUSIC CALC-BMR
#-----------------------------------------------------------------------

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

dia-high-af/gene_mrs: dia-high-af/total_covgs dia/wig-list.tsv dia-high-af/rem_dia.maf ~/hdall/data/hg19/ucsc.hg19.nochr.fasta ucsc-genes.hg19.roi
	rm -f dia-high-af/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list dia/wig-list.tsv \
		--maf-file dia-high-af/rem_dia.maf \
		--output-dir dia-high-af/ \
		--reference-sequence ~/hdall/data/hg19/ucsc.hg19.nochr.fasta \
		--roi-file ucsc-genes.hg19.roi

rel-high-af/gene_mrs: rel-high-af/total_covgs rel/wig-list.tsv rel-high-af/rem_rel.maf ~/hdall/data/hg19/ucsc.hg19.nochr.fasta ucsc-genes.hg19.roi
	rm -f rel-high-af/overall_bmrs
	genome music bmr calc-bmr \
		--bam-list rel/wig-list.tsv \
		--maf-file rel-high-af/rem_rel.maf \
		--output-dir rel-high-af/ \
		--reference-sequence ~/hdall/data/hg19/ucsc.hg19.nochr.fasta \
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

dia-high-af/smg.tsv: dia-high-af/gene_mrs
	genome music smg \
		--gene-mr-file dia-high-af/gene_mrs \
		--output-file dia-high-af/smg.tsv.part
	mv dia-high-af/smg.tsv.part dia-high-af/smg.tsv

rel-high-af/smg.tsv: rel-high-af/gene_mrs
	genome music smg \
		--gene-mr-file rel-high-af/gene_mrs \
		--output-file rel-high-af/smg.tsv.part
	mv rel-high-af/smg.tsv.part rel-high-af/smg.tsv

#-----------------------------------------------------------------------
#--- GENOME MUSIC PATH-SCAN (using only high allelic frequency variants)
#-----------------------------------------------------------------------

pathscan/pathways-david-complete.tsv: pathscan/pathways-david.tsv pathscan/pathways-david-reverse.tsv
	perl ~/hdall/scripts/merge-david-annotations.pl \
		pathscan/pathways-david.tsv \
		pathscan/pathways-david-reverse.tsv \
		> pathscan/pathways-david-complete.tsv.part
	mv pathscan/pathways-david-complete.tsv.part pathscan/pathways-david-complete.tsv

pathscan/gene-sets.mutsigdb4.0.tsv: ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt ~/hdall/scripts/msigdb2pathscan.pl
	cat ~/hdall/data/msigdb/c1.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c1_positional \
		2>&1 1>pathscan/gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c2.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c2_curated \
		2>&1 1>>pathscan/gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c3.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c3_motif \
		2>&1 1>>pathscan/gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c4.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c4_computational \
		2>&1 1>>pathscan/gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c5.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c5_GO \
		2>&1 1>>pathscan/gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c6.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c6_oncogenic_signature \
		2>&1 1>>pathscan/gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c7.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c7_immunologic_signature \
		2>&1 1>>pathscan/gene-sets.mutsigdb4.0.tsv.part | tee -a make.log

pathscan/pathways-nci.tsv: ~/hdall/data/nci-pathway/uniprot.tab ~/hdall/scripts/nci2pathscan.pl
	cat ~/hdall/data/nci-pathway/uniprot.tab | perl ~/hdall/scripts/nci2pathscan.pl \
		2>&1 1>>pathscan/pathways-nci.tsv.part | tee -a make.log
	mv pathscan/pathways-nci.tsv.part pathscan/pathways-nci.tsv

pathscan/pathways-wikipathways.tsv: ~/hdall/data/wikipathways/wikipathways_Homo_sapiens_Curation-AnalysisCollection__txt.zip
	perl ~/hdall/scripts/wikipathways2pathscan.pl \
		--wiki-dir ~/hdall/data/wikipathways/ \
		2>&1 1>>pathscan/pathways-wikipathways.tsv.part | tee -a make.log
	mv pathscan/pathways-wikipathways.tsv.part pathscan/pathways-wikipathways.tsv
		
pathscan/pathways-all.tsv: pathscan/pathways-david-complete.tsv pathscan/gene-sets.mutsigdb4.0.tsv pathscan/pathways-nci.tsv pathscan/pathways-wikipathways.tsv
	cp pathscan/pathways-david-complete.tsv pathscan/pathways-all.tsv.part
	cat pathscan/gene-sets.mutsigdb4.0.tsv pathscan/pathways-nci.tsv pathscan/pathways-manual.tsv pathscan/pathways-wikipathways.tsv >> pathscan/pathways-all.tsv.part
	mv pathscan/pathways-all.tsv.part pathscan/pathways-all.tsv

pathscan/pathways-all.sizefiltered.tsv: pathscan/pathways-all.tsv
	cat pathscan/pathways-all.tsv | perl -ne '@f = split /\t/; @g = split(/\|/, $$f[3]); print $$_ if (@g > 1 and @g < 400);' > pathscan/pathways-all.sizefiltered.tsv.part
	mv pathscan/pathways-all.sizefiltered.tsv.part pathscan/pathways-all.sizefiltered.tsv

pathscan/pathways-kegg-biocarta-bbid-reactome.tsv: pathscan/pathways-david-complete.tsv pathscan/gene-sets.mutsigdb4.0.tsv
	grep -P "(KEGG_PATHWAY|BIOCARTA|BBID)" pathscan/pathways-david-complete.tsv > pathscan/pathways-kegg-biocarta-bbid-reactome.tsv
	grep "REACTOME_" pathscan/gene-sets.mutsigdb4.0.tsv >> pathscan/pathways-kegg-biocarta-bbid-reactome.tsv


# iteration 1: exclude likely false positives from Lawrence et al. (2013)
#-------------------------------------------------------------------------
pathscan/sm_pathways.dia.tsv: dia-high-af/smg.tsv pathscan/pathways-all.sizefiltered.tsv
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia-high-af/gene_covgs/ \
		--maf-file dia-high-af/rem_dia.maf \
		--output-file pathscan/sm_pathways.dia.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1
	mv pathscan/sm_pathways.dia.tsv.part pathscan/sm_pathways.dia.tsv
	mv pathscan/sm_pathways.dia.tsv.part_detailed pathscan/sm_pathways.dia.tsv_detailed

pathscan/sm_pathways.dia.annotated.tsv: pathscan/sm_pathways.dia.tsv ~/hdall/scripts/music/annotate-pathscan-result.pl
	perl ~/hdall/scripts/music/annotate-pathscan-result.pl \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--sm-pathways pathscan/sm_pathways.dia.tsv \
		--sm-pathways-detail pathscan/sm_pathways.dia.tsv_detailed \
		--gene-pathway-matrix pathscan/sm_pathways.dia.matrix.tsv \
		> pathscan/sm_pathways.dia.annotated.tsv.part
	mv pathscan/sm_pathways.dia.annotated.tsv.part pathscan/sm_pathways.dia.annotated.tsv

pathscan/sm_pathways.rel.tsv: rel-high-af/smg.tsv pathscan/pathways-all.sizefiltered.tsv
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel-high-af/gene_covgs/ \
		--maf-file rel-high-af/rem_rel.maf \
		--output-file pathscan/sm_pathways.rel.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1
	mv pathscan/sm_pathways.rel.tsv.part pathscan/sm_pathways.rel.tsv
	mv pathscan/sm_pathways.rel.tsv.part_detailed pathscan/sm_pathways.rel.tsv_detailed

pathscan/sm_pathways.rel.annotated.tsv: pathscan/sm_pathways.rel.tsv ~/hdall/scripts/music/annotate-pathscan-result.pl
	perl ~/hdall/scripts/music/annotate-pathscan-result.pl \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--sm-pathways pathscan/sm_pathways.rel.tsv \
		--sm-pathways-detail pathscan/sm_pathways.rel.tsv_detailed \
		--gene-pathway-matrix pathscan/sm_pathways.rel.matrix.tsv \
		> pathscan/sm_pathways.rel.annotated.tsv.part
	mv pathscan/sm_pathways.rel.annotated.tsv.part pathscan/sm_pathways.rel.annotated.tsv

pathscan/sm_pathways.tsv: pathscan/sm_pathways.dia.annotated.tsv pathscan/sm_pathways.rel.annotated.tsv ~/hdall/scripts/cluster-pathways.R ~/hdall/scripts/music/merge-sm-pathways.R
	R --no-save --quiet --slave -f ~/hdall/scripts/music/merge-sm-pathways.R \
		--args pathscan/sm_pathways.dia.annotated.tsv pathscan/sm_pathways.rel.annotated.tsv \
		2>&1 1>pathscan/sm_pathways.tsv.part | tee -a make.log
	R --no-save --quiet --slave -f ~/hdall/scripts/cluster-pathways.R \
		--args pathscan/sm_pathways.tsv.part pathscan/sm_pathways.clustered.tsv.part 13 1e-7 \
		2>&1 | tee -a make.log
	mv pathscan/sm_pathways.clustered.tsv.part pathscan/sm_pathways.tsv
	rm pathscan/sm_pathways.tsv.part

../pathway-patient-matrix.clustered.tsv: pathscan/sm_pathways.dia.annotated.tsv pathscan/sm_pathways.rel.annotated.tsv ../filtered-variants.tsv ../kamilla/gene-pathway-assignment.tsv ~/hdall/scripts/get-pathway-patient-matrix.pl ~/hdall/scripts/cluster-pathways.R
	perl ~/hdall/scripts/get-pathway-patient-matrix.pl \
		--enriched-pathways-dia pathscan/sm_pathways.dia.annotated.tsv \
		--enriched-pathways-rel pathscan/sm_pathways.rel.annotated.tsv \
		--curated-pathways ../kamilla/gene-pathway-assignment.tsv \
		--filtered-variants ../filtered-variants.tsv \
		2>&1 1>../pathway-patient-matrix.tsv.part | tee -a make.log
	mv ../pathway-patient-matrix.tsv.part ../pathway-patient-matrix.tsv
	# TABLE: pathway-patient-matrix
	R --no-save --quiet --slave -f ~/hdall/scripts/cluster-pathways.R \
		--args ../pathway-patient-matrix.tsv ../pathway-patient-matrix.clustered.tsv.part 27 1e-7 \
		2>&1 | tee -a make.log
	mv ../pathway-patient-matrix.clustered.tsv.part ../pathway-patient-matrix.clustered.tsv

# iteration 2: additionally exclude KRAS, NRAS, and CREBBP
#----------------------------------------------------------
pathscan/sm_pathways.dia.iter2.tsv: dia-high-af/smg.tsv pathscan/pathways-all.sizefiltered.tsv
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia-high-af/gene_covgs/ \
		--maf-file dia-high-af/rem_dia.maf \
		--output-file pathscan/sm_pathways.dia.iter2.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore KRAS,NRAS,CREBBP,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1
	mv pathscan/sm_pathways.dia.iter2.tsv.part pathscan/sm_pathways.dia.iter2.tsv
	mv pathscan/sm_pathways.dia.iter2.tsv.part_detailed pathscan/sm_pathways.dia.iter2.tsv_detailed

pathscan/sm_pathways.dia.iter2.annotated.tsv: pathscan/sm_pathways.dia.iter2.tsv ~/hdall/scripts/music/annotate-pathscan-result.pl
	perl ~/hdall/scripts/music/annotate-pathscan-result.pl \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--sm-pathways pathscan/sm_pathways.dia.iter2.tsv \
		--sm-pathways-detail pathscan/sm_pathways.dia.iter2.tsv_detailed \
		--gene-pathway-matrix pathscan/sm_pathways.dia.iter2.matrix.tsv \
		> pathscan/sm_pathways.dia.iter2.annotated.tsv.part
	mv pathscan/sm_pathways.dia.iter2.annotated.tsv.part pathscan/sm_pathways.dia.iter2.annotated.tsv

pathscan/sm_pathways.rel.iter2.tsv: rel-high-af/smg.tsv pathscan/pathways-all.sizefiltered.tsv
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel-high-af/gene_covgs/ \
		--maf-file rel-high-af/rem_rel.maf \
		--output-file pathscan/sm_pathways.rel.iter2.tsv.part \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--genes-to-ignore KRAS,NRAS,CREBBP,TTN,TBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1
	mv pathscan/sm_pathways.rel.iter2.tsv.part pathscan/sm_pathways.rel.iter2.tsv
	mv pathscan/sm_pathways.rel.iter2.tsv.part_detailed pathscan/sm_pathways.rel.iter2.tsv_detailed

pathscan/sm_pathways.rel.iter2.annotated.tsv: pathscan/sm_pathways.rel.iter2.tsv ~/hdall/scripts/music/annotate-pathscan-result.pl
	perl ~/hdall/scripts/music/annotate-pathscan-result.pl \
		--pathway-file pathscan/pathways-all.sizefiltered.tsv \
		--sm-pathways pathscan/sm_pathways.rel.iter2.tsv \
		--sm-pathways-detail pathscan/sm_pathways.rel.iter2.tsv_detailed \
		--gene-pathway-matrix pathscan/sm_pathways.rel.iter2.matrix.tsv \
		> pathscan/sm_pathways.rel.iter2.annotated.tsv.part
	mv pathscan/sm_pathways.rel.iter2.annotated.tsv.part pathscan/sm_pathways.rel.iter2.annotated.tsv

pathscan/sm_pathways.iter2.tsv: pathscan/sm_pathways.dia.iter2.annotated.tsv pathscan/sm_pathways.rel.iter2.annotated.tsv ~/hdall/scripts/cluster-pathways.R ~/hdall/scripts/music/merge-sm-pathways.R
	R --no-save --quiet --slave -f ~/hdall/scripts/music/merge-sm-pathways.R \
		--args pathscan/sm_pathways.dia.iter2.annotated.tsv pathscan/sm_pathways.rel.iter2.annotated.tsv \
		2>&1 1>pathscan/sm_pathways.iter2.tsv.part | tee -a make.log
	R --no-save --quiet --slave -f ~/hdall/scripts/cluster-pathways.R \
		--args pathscan/sm_pathways.iter2.tsv.part pathscan/sm_pathways.iter2.clustered.tsv.part 13 1e-7 \
		2>&1 | tee -a make.log
	mv pathscan/sm_pathways.iter2.clustered.tsv.part pathscan/sm_pathways.iter2.tsv
	rm pathscan/sm_pathways.iter2.tsv.part

../pathway-patient-matrix.iter2.clustered.tsv: pathscan/sm_pathways.dia.iter2.annotated.tsv pathscan/sm_pathways.rel.iter2.annotated.tsv pathscan/sm_pathways.iter2.tsv ../filtered-variants.tsv ../kamilla/gene-pathway-assignment.tsv ~/hdall/scripts/get-pathway-patient-matrix.pl ~/hdall/scripts/cluster-pathways.R
	perl ~/hdall/scripts/get-pathway-patient-matrix.pl \
		--enriched-pathways-dia pathscan/sm_pathways.dia.iter2.annotated.tsv \
		--enriched-pathways-rel pathscan/sm_pathways.rel.iter2.annotated.tsv \
		--curated-pathways ../kamilla/gene-pathway-assignment.tsv \
		--filtered-variants ../filtered-variants.tsv \
		2>&1 1>../pathway-patient-matrix.iter2.tsv.part | tee -a make.log
	mv ../pathway-patient-matrix.iter2.tsv.part ../pathway-patient-matrix.iter2.tsv
	# TABLE: pathway-patient-matrix
	R --no-save --quiet --slave -f ~/hdall/scripts/cluster-pathways.R \
		--args ../pathway-patient-matrix.iter2.tsv ../pathway-patient-matrix.iter2.clustered.tsv.part 27 1e-7 \
		2>&1 | tee -a make.log
	mv ../pathway-patient-matrix.iter2.clustered.tsv.part ../pathway-patient-matrix.iter2.clustered.tsv

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

pathscan/gene-pathway-matrix.clustered.tsv: pathscan/pathways-nci.tsv pathscan/pathways-kegg-biocarta-bbid-reactome.tsv pathscan/pathways-wikipathways.tsv ~/hdall/scripts/get-gene-pathway-matrix.pl ~/hdall/scripts/cluster-gene-pathway-matrix.R 
	cat pathscan/pathways-nci.tsv pathscan/pathways-kegg-biocarta-bbid-reactome.tsv pathscan/pathways-wikipathways.tsv \
		| perl ~/hdall/scripts/get-gene-pathway-matrix.pl --pathway-file - > pathscan/gene-pathway-matrix.tsv.part
	mv pathscan/gene-pathway-matrix.tsv.part pathscan/gene-pathway-matrix.tsv
	R --no-save --quiet --slave -f ~/hdall/scripts/cluster-gene-pathway-matrix.R \
		--args pathscan/gene-pathway-matrix.tsv pathscan/gene-pathway-matrix.clustered.tsv \
		2>&1 | tee -a make.log

../gene-patient-matrix.annotated.tsv: ../gene-patient-matrix.tsv dia/smg.tsv rel/smg.tsv pathscan/sm_pathways.dia.annotated.tsv pathscan/sm_pathways.rel.annotated.tsv ../cnv/impacted-genes-list.cnv.tsv pathscan/gene-pathway-matrix.clustered.tsv ~/hdall/scripts/annotate-gene-patient-matrix.pl
	perl ~/hdall/scripts/annotate-gene-patient-matrix.pl \
		--gene-patient-matrix ../gene-patient-matrix.tsv \
		--smg-dia dia/smg.tsv \
		--smg-rel rel/smg.tsv \
		--smp-dia pathscan/sm_pathways.dia.annotated.tsv \
		--smp-rel pathscan/sm_pathways.rel.annotated.tsv \
		--cnv ../cnv/impacted-genes-list.cnv.tsv \
		--gene-pathway-matrix pathscan/gene-pathway-matrix.clustered.tsv \
		2>&1 1>../gene-patient-matrix.annotated.tsv.part | tee -a make.log
	mv ../gene-patient-matrix.annotated.tsv.part ../gene-patient-matrix.annotated.tsv	
	
../gene-patient-matrix.high-af.annotated.tsv: ../gene-patient-matrix.high-af.tsv dia-high-af/smg.tsv rel-high-af/smg.tsv pathscan/sm_pathways.dia.annotated.tsv pathscan/sm_pathways.rel.annotated.tsv ../cnv/impacted-genes-list.cnv.tsv pathscan/gene-pathway-matrix.clustered.tsv ~/hdall/scripts/annotate-gene-patient-matrix.pl
	perl ~/hdall/scripts/annotate-gene-patient-matrix.pl \
		--gene-patient-matrix ../gene-patient-matrix.high-af.tsv \
		--smg-dia dia-high-af/smg.tsv \
		--smg-rel rel-high-af/smg.tsv \
		--smp-dia pathscan/sm_pathways.dia.annotated.tsv \
		--smp-rel pathscan/sm_pathways.rel.annotated.tsv \
		--cnv ../cnv/impacted-genes-list.cnv.tsv \
		--gene-pathway-matrix pathscan/gene-pathway-matrix.clustered.tsv \
		2>&1 1>../gene-patient-matrix.high-af.annotated.tsv.part | tee -a make.log
	mv ../gene-patient-matrix.high-af.annotated.tsv.part ../gene-patient-matrix.high-af.annotated.tsv
		
	