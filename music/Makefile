export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work

all: dia rel ../gene-patient-matrix.annotated.tsv ../gene-patient-matrix.high-af.annotated.tsv sm_pathways.tsv sm_pathways.high-af.tsv MSigDB rel/mutation-relation.tsv ../pathway-patient-matrix.tsv ../pathway-patient-matrix.high-af.tsv

dia: ucsc-genes.hg19.roi dia/rem_dia.maf dia/gene_mrs dia/smg.tsv dia/sm_pathways.tsv dia/sm_pathways.annotated.tsv
rel: ucsc-genes.hg19.roi rel/rem_rel.maf rel/gene_mrs rel/smg.tsv rel/sm_pathways.tsv rel/sm_pathways.annotated.tsv
dia-high-af: ucsc-genes.hg19.roi dia-high-af/rem_dia.maf dia-high-af/gene_mrs dia-high-af/smg.tsv dia-high-af/sm_pathways.tsv dia-high-af/sm_pathways.annotated.tsv
rel-high-af: ucsc-genes.hg19.roi rel-high-af/rem_rel.maf rel-high-af/gene_mrs rel-high-af/smg.tsv rel-high-af/sm_pathways.tsv rel-high-af/sm_pathways.annotated.tsv

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
		2>&1 1>$@.part | grep -vP '(Leading or trailing space|variant.Format)' | tee -a make.log
	cat ../filtered_variants/$*.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		2>&1 1>>$@.part | grep -vP '(Leading or trailing space|variant.Format)' | tee -a make.log
	mv $@.part $@

maf-high-af/%.maf: ../filtered_variants/%.snp.filtered.vcf ../filtered_variants/%.indel.filtered.vcf ucsc-genes.hg19.roi ~/hdall/results/gene-id-mapping.biomart-0.7.tsv ~/hdall/scripts/vcf2maf.pl
	mkdir -p maf-high-af
	cat ../filtered_variants/$*.snp.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.25 \
		2>&1 1>$@.part | grep -vP '(Leading or trailing space|variant.Format)' | tee -a make.log
	cat ../filtered_variants/$*.indel.filtered.vcf | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi ucsc-genes.hg19.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.25 \
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

#--- GENOME MUSIC PROXIMITY

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
	mv dia/sm_pathways.tsv.part_detailed dia/sm_pathways.tsv_detailed

dia/sm_pathways.annotated.tsv: dia/sm_pathways.tsv
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
	mv rel/sm_pathways.tsv.part_detailed rel/sm_pathways.tsv_detailed

rel/sm_pathways.annotated.tsv: rel/sm_pathways.tsv
	perl ~/hdall/scripts/music/annotate-pathscan-result.pl \
		--sm-pathways rel/sm_pathways.tsv \
		--sm-pathways-detail rel/sm_pathways.tsv_detailed \
		> rel/sm_pathways.annotated.tsv.part
	mv rel/sm_pathways.annotated.tsv.part rel/sm_pathways.annotated.tsv

dia-high-af/sm_pathways.tsv: dia-high-af/smg.tsv pathways-david-complete.tsv
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia-high-af/gene_covgs/ \
		--maf-file dia-high-af/rem_dia.maf \
		--genes-to-ignore TTN,PTPN11,NF1,TBP,KRAS,NRAS,CREBBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1 \
		--output-file dia-high-af/sm_pathways.tsv.part \
		--pathway-file pathways-david-complete.tsv
	mv dia-high-af/sm_pathways.tsv.part dia-high-af/sm_pathways.tsv
	mv dia-high-af/sm_pathways.tsv.part_detailed dia-high-af/sm_pathways.tsv_detailed

dia-high-af/sm_pathways.annotated.tsv: dia-high-af/sm_pathways.tsv
	perl ~/hdall/scripts/music/annotate-pathscan-result.pl \
		--sm-pathways dia-high-af/sm_pathways.tsv \
		--sm-pathways-detail dia-high-af/sm_pathways.tsv_detailed \
		> dia-high-af/sm_pathways.annotated.tsv.part
	mv dia-high-af/sm_pathways.annotated.tsv.part dia-high-af/sm_pathways.annotated.tsv

rel-high-af/sm_pathways.tsv: rel-high-af/smg.tsv pathways-david-complete.tsv
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel-high-af/gene_covgs/ \
		--maf-file rel-high-af/rem_rel.maf \
		--genes-to-ignore TTN,PTPN11,NF1,TBP,KRAS,NRAS,CREBBP,CSMD3,DNAH5,RYR1,RYR2,RYR3,DNAH1,DNAH8,DNAH9,MUC2,MUC16,MUC12,MUC5B,OR5H2,OR5H2,OR11H4,OR6F1,OR52R1,OR2T12,OR6V1,OR51I2,OR5I1,OR9Q1,OR9Q1 \
		--output-file rel-high-af/sm_pathways.tsv.part \
		--pathway-file pathways-david-complete.tsv
	mv rel-high-af/sm_pathways.tsv.part rel-high-af/sm_pathways.tsv
	mv rel-high-af/sm_pathways.tsv.part_detailed rel-high-af/sm_pathways.tsv_detailed

rel-high-af/sm_pathways.annotated.tsv: rel-high-af/sm_pathways.tsv
	perl ~/hdall/scripts/music/annotate-pathscan-result.pl \
		--sm-pathways rel-high-af/sm_pathways.tsv \
		--sm-pathways-detail rel-high-af/sm_pathways.tsv_detailed \
		> rel-high-af/sm_pathways.annotated.tsv.part
	mv rel-high-af/sm_pathways.annotated.tsv.part rel-high-af/sm_pathways.annotated.tsv

sm_pathways.tsv: dia/sm_pathways.annotated.tsv rel/sm_pathways.annotated.tsv ~/hdall/scripts/cluster-pathways.R ~/hdall/scripts/music/merge-sm-pathways.R
	R --no-save --quiet --slave -f ~/hdall/scripts/music/merge-sm-pathways.R \
		2>&1 1>sm_pathways.tsv.part | tee -a make.log
	R --no-save --quiet --slave -f ~/hdall/scripts/cluster-pathways.R \
		--args sm_pathways.tsv.part sm_pathways.clustered.tsv.part \
		2>&1 | tee -a make.log
	mv sm_pathways.clustered.tsv.part sm_pathways.tsv
	rm sm_pathways.tsv.part

sm_pathways.high-af.tsv: dia-high-af/sm_pathways.annotated.tsv rel-high-af/sm_pathways.annotated.tsv ~/hdall/scripts/cluster-pathways.R ~/hdall/scripts/music/merge-sm-pathways.R
	R --no-save --quiet --slave -f ~/hdall/scripts/music/merge-sm-pathways.R \
		2>&1 1>sm_pathways.high-af.tsv.part | tee -a make.log
	R --no-save --quiet --slave -f ~/hdall/scripts/cluster-pathways.R \
		--args sm_pathways.high-af.tsv.part sm_pathways.high-af.clustered.tsv.part \
		2>&1 | tee -a make.log
	mv sm_pathways.high-af.clustered.tsv.part sm_pathways.high-af.tsv
	rm sm_pathways.high-af.tsv.part

#--- MSigDB gene set enrichment
	
MSigDB: gene-sets.mutsigdb4.0.tsv dia/sm_pathways.msigdb.annotated.tsv rel/sm_pathways.msigdb.annotated.tsv

gene-sets.mutsigdb4.0.tsv: ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt ~/hdall/scripts/msigdb2pathscan.pl
	cat ~/hdall/data/msigdb/c1.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c1_positional \
		2>&1 1>gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c2.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c2_curated \
		2>&1 1>>gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c3.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c3_motif \
		2>&1 1>>gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c4.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c4_computational \
		2>&1 1>>gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c5.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c5_GO \
		2>&1 1>>gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c6.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c6_oncogenic_signature \
		2>&1 1>>gene-sets.mutsigdb4.0.tsv.part | tee -a make.log 
	cat ~/hdall/data/msigdb/c7.all.v4.0.symbols.gmt | perl ~/hdall/scripts/msigdb2pathscan.pl \
		--entrez-id-file ~/hdall/data/msigdb/msigdb.v4.0.entrez.gmt \
		--category c7_immunologic_signature \
		2>&1 1>>gene-sets.mutsigdb4.0.tsv.part | tee -a make.log
	mv gene-sets.mutsigdb4.0.tsv.part gene-sets.mutsigdb4.0.tsv

dia/sm_pathways.msigdb.tsv: dia/smg.tsv dia/wig-list.tsv dia-high-af/rem_dia.maf gene-sets.mutsigdb4.0.tsv
	genome music path-scan \
		--bam-list dia/wig-list.tsv \
		--gene-covg-dir dia-high-af/gene_covgs/ \
		--maf-file dia-high-af/rem_dia.maf \
		--output-file dia/sm_pathways.msigdb.tsv.part \
		--pathway-file gene-sets.mutsigdb4.0.tsv
	#	--noskip-non-coding \
	#	--noskip-silent \
	mv dia/sm_pathways.msigdb.tsv.part dia/sm_pathways.msigdb.tsv
	mv dia/sm_pathways.msigdb.tsv.part_detailed dia/sm_pathways.msigdb.tsv_detailed

dia/sm_pathways.msigdb.annotated.tsv: dia/sm_pathways.msigdb.tsv
	perl ~/hdall/scripts/music/annotate-pathscan-result.pl \
		--sm-pathways dia/sm_pathways.msigdb.tsv \
		--sm-pathways-detail dia/sm_pathways.msigdb.tsv_detailed \
		> dia/sm_pathways.msigdb.annotated.tsv.part
	mv dia/sm_pathways.msigdb.annotated.tsv.part dia/sm_pathways.msigdb.annotated.tsv

rel/sm_pathways.msigdb.tsv: rel/smg.tsv rel/wig-list.tsv rel-high-af/rem_rel.maf gene-sets.mutsigdb4.0.tsv
	genome music path-scan \
		--bam-list rel/wig-list.tsv \
		--gene-covg-dir rel-high-af/gene_covgs/ \
		--maf-file rel-high-af/rem_rel.maf \
		--output-file rel/sm_pathways.msigdb.tsv.part \
		--pathway-file gene-sets.mutsigdb4.0.tsv
	#	--noskip-non-coding \
	#	--noskip-silent \
	mv rel/sm_pathways.msigdb.tsv.part rel/sm_pathways.msigdb.tsv
	mv rel/sm_pathways.msigdb.tsv.part_detailed rel/sm_pathways.msigdb.tsv_detailed

rel/sm_pathways.msigdb.annotated.tsv: rel/sm_pathways.msigdb.tsv
	perl ~/hdall/scripts/music/annotate-pathscan-result.pl \
		--sm-pathways rel/sm_pathways.msigdb.tsv \
		--sm-pathways-detail rel/sm_pathways.msigdb.tsv_detailed \
		> rel/sm_pathways.msigdb.annotated.tsv.part
	mv rel/sm_pathways.msigdb.annotated.tsv.part rel/sm_pathways.msigdb.annotated.tsv

#--- MUSIC MUTATION-RELATION ANALYSIS

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

#--- COMBINE ALL RESULTS INTO FINAL LIST

../gene-patient-matrix.annotated.tsv: ../gene-patient-matrix.tsv dia/smg.tsv rel/smg.tsv dia/sm_pathways.annotated.tsv rel/sm_pathways.annotated.tsv ~/hdall/scripts/annotate-gene-patient-matrix.pl
	perl ~/hdall/scripts/annotate-gene-patient-matrix.pl \
		--gene-patient-matrix ../gene-patient-matrix.tsv \
		--smg-dia dia/smg.tsv \
		--smg-rel rel/smg.tsv \
		--smp-dia dia/sm_pathways.annotated.tsv \
		--smp-rel rel/sm_pathways.annotated.tsv \
		--cnv ../cnv/impacted-genes-list.cnv.tsv \
		2>&1 1>../gene-patient-matrix.annotated.tsv.part | tee -a make.log
	mv ../gene-patient-matrix.annotated.tsv.part ../gene-patient-matrix.annotated.tsv	
	
../gene-patient-matrix.high-af.annotated.tsv: ../gene-patient-matrix.high-af.tsv dia-high-af/smg.tsv rel-high-af/smg.tsv dia-high-af/sm_pathways.annotated.tsv rel-high-af/sm_pathways.annotated.tsv ~/hdall/scripts/annotate-gene-patient-matrix.pl
	perl ~/hdall/scripts/annotate-gene-patient-matrix.pl \
		--gene-patient-matrix ../gene-patient-matrix.high-af.tsv \
		--smg-dia dia-high-af/smg.tsv \
		--smg-rel rel-high-af/smg.tsv \
		--smp-dia dia-high-af/sm_pathways.annotated.tsv \
		--smp-rel rel-high-af/sm_pathways.annotated.tsv \
		--cnv ../cnv/impacted-genes-list.cnv.tsv \
		2>&1 1>../gene-patient-matrix.high-af.annotated.tsv.part | tee -a make.log
	mv ../gene-patient-matrix.high-af.annotated.tsv.part ../gene-patient-matrix.high-af.annotated.tsv
	
../pathway-patient-matrix.tsv: dia/sm_pathways.annotated.tsv rel/sm_pathways.annotated.tsv ../gene-patient-matrix.tsv ~/hdall/scripts/get-pathway-patient-matrix.pl ~/hdall/scripts/cluster-pathways.R
	perl ~/hdall/scripts/get-pathway-patient-matrix.pl \
		--enriched-pathways-dia dia/sm_pathways.annotated.tsv \
		--enriched-pathways-rel rel/sm_pathways.annotated.tsv \
		--gene-patient-matrix ../gene-patient-matrix.tsv \
		2>&1 1>../pathway-patient-matrix.tsv.part | tee -a make.log
	R --no-save --quiet --slave -f ~/hdall/scripts/cluster-pathways.R \
		--args ../pathway-patient-matrix.tsv.part ../pathway-patient-matrix.clustered.tsv
	mv ../pathway-patient-matrix.clustered.tsv ../pathway-patient-matrix.tsv
	rm ../pathway-patient-matrix.tsv.part

../pathway-patient-matrix.high-af.tsv: dia-high-af/sm_pathways.annotated.tsv rel-high-af/sm_pathways.annotated.tsv ../gene-patient-matrix.high-af.tsv ~/hdall/scripts/get-pathway-patient-matrix.pl ~/hdall/scripts/cluster-pathways.R
	perl ~/hdall/scripts/get-pathway-patient-matrix.pl \
		--enriched-pathways-dia dia-high-af/sm_pathways.annotated.tsv \
		--enriched-pathways-rel rel-high-af/sm_pathways.annotated.tsv \
		--gene-patient-matrix ../gene-patient-matrix.high-af.tsv \
		2>&1 1>../pathway-patient-matrix.high-af.tsv.part | tee -a make.log
	R --no-save --quiet --slave -f ~/hdall/scripts/cluster-pathways.R \
		--args ../pathway-patient-matrix.high-af.tsv.part ../pathway-patient-matrix.high-af.clustered.tsv
	mv ../pathway-patient-matrix.high-af.clustered.tsv ../pathway-patient-matrix.high-af.tsv
	rm ../pathway-patient-matrix.high-af.tsv.part
	
	