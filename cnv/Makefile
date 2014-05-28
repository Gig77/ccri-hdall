export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

PATIENTS_EXOME = 314 399 430 446 460 545 592 715 786 792 818 842 A B C D X Y
NOT_WORKING = 1021247 E
PATIENTS_RESEQ_MATCHED = 1009302 1019357 1020540 1020583 1021392 1021631 1022914 1023056 1023392 1024518 1024543 1025108 1025409 1187 FB11 G HS6 K MB2 243 933 944 KA14651 KD20493 KE12025 MJ16441 NH17331 PJ13414 RD17412 RS13466 ST13892 ZA16211
PATIENTS_RESEQ_RELONLY = 1017005 1021865 1023545 AD15 BL16 BM18 CA18 DM1 FE1 FS1 GD18 GD1 HJ15 HJA15 HL1 KA17 KJ17 KL16 LB17 LM18 MJ1 ML10 MV16 NS18 PC16 RT15 RT16 SJM16 SKR1 SL1 SLM1 ST14 WA1 ZE13
PATIENTS_RESEQ_DIAONLY = 1004564 1010661 1010781 1019964 1020076 1021087 1023338 1023616 1024589 1026233 1026662 B100 EF7 FB14 G44 HD7
PATIENTS_RESEQ_NONREL = 331 380 442 350 461 466 529 591 602 619 633 634 642 653 666 672 697 698 700 709 724 762 776 777 779 782 409 NRD_1 73 NRD_2 NRD_3 60 594 687 748 754 646 530 718 681 39 49 45 54 110 111 134 143 147 199 NRD_4
#PATIENTS_EXOME = 
#PATIENTS_RESEQ_MATCHED =  
#PATIENTS_RESEQ_RELONLY = 
#PATIENTS_RESEQ_DIAONLY = 


all: $(foreach P, $(PATIENTS_EXOME), $P_rem.coverage.exome.tsv $P_dia.coverage.exome.tsv $P_rel.coverage.exome.tsv $P_rem.coverage.panel.tsv $P_dia.coverage.panel.tsv $P_rel.coverage.panel.tsv ../reseq/cnv/$P_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/$P_Relapse.coverage.reseq.panel.tsv cov-plot.$P.exome.pdf cov-plot.$P.panel.pdf gc/$P_rem.coverage.bed gc/$P_dia.coverage.bed gc/$P_rel.coverage.bed) \
	 $(foreach P, $(PATIENTS_RESEQ_MATCHED), ../reseq/cnv/$P_Remission.coverage.reseq.panel.tsv ../reseq/cnv/$P_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/$P_Relapse.coverage.reseq.panel.tsv ../reseq/cnv/cov-plot.$P.matched.reseq.panel.pdf ../reseq/cnv/segmented_coverage/$P_Diagnosis.segmented-coverage.tsv  ../reseq/cnv/segmented_coverage/$P_Remission.segmented-coverage.tsv  ../reseq/cnv/segmented_coverage/$P_Relapse.segmented-coverage.tsv ../reseq/cnv/segmented_coverage/$P_Diagnosis.segmented-coverage.pdf ../reseq/cnv/segmented_coverage/$P_Relapse.segmented-coverage.pdf) \
	 $(foreach P, $(PATIENTS_RESEQ_RELONLY), ../reseq/cnv/$P_Remission.coverage.reseq.panel.tsv ../reseq/cnv/$P_Relapse.coverage.reseq.panel.tsv ../reseq/cnv/cov-plot.$P.rel.reseq.panel.pdf ../reseq/cnv/segmented_coverage/$P_Remission.segmented-coverage.tsv  ../reseq/cnv/segmented_coverage/$P_Relapse.segmented-coverage.tsv ../reseq/cnv/segmented_coverage/$P_Relapse.segmented-coverage.pdf) \
	 $(foreach P, $(PATIENTS_RESEQ_DIAONLY), ../reseq/cnv/$P_Remission.coverage.reseq.panel.tsv ../reseq/cnv/$P_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/cov-plot.$P.dia.reseq.panel.pdf ../reseq/cnv/segmented_coverage/$P_Diagnosis.segmented-coverage.tsv  ../reseq/cnv/segmented_coverage/$P_Remission.segmented-coverage.tsv ../reseq/cnv/segmented_coverage/$P_Diagnosis.segmented-coverage.pdf) \
	 $(foreach P, $(PATIENTS_RESEQ_NONREL), ../reseq/cnv/$P_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/cov-plot.$P.dia.nonrel.reseq.panel.pdf ../reseq/cnv/segmented_coverage/$P_Diagnosis.segmented-coverage.tsv ../reseq/cnv/segmented_coverage/$P_Diagnosis.nonrel.segmented-coverage.pdf) \
	 relapsing.rel.circos.png

# get GC content for each exon
gc/truseq_exome_targeted_regions.hg19.gc.bed: ~/generic/data/hg19/ucsc.hg19.fasta ~/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed.chr
	bedtools nuc \
		-fi ~/generic/data/hg19/ucsc.hg19.fasta \
		-bed ~/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed.chr \
		| cut -f 1-4,8 > $@.part
	mv $@.part $@

# append read depth
gc/%.coverage.bed: ~/hdall/data/bam/%.merged.duplicate_marked.realigned.recalibrated.bam gc/truseq_exome_targeted_regions.hg19.gc.bed
	bedtools coverage \
		-abam ~/hdall/data/bam/$*.merged.duplicate_marked.realigned.recalibrated.bam \
		-b gc/truseq_exome_targeted_regions.hg19.gc.bed \
		> $@.part
	mv $@.part $@

%.coverage.panel.tsv: /data/christian/hdall/data/current/bam/%.merged.duplicate_marked.realigned.recalibrated.bam /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr 
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b ~/hdall/data/reseq/enriched_regions.bed \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed ~/hdall/data/reseq/enriched_regions.bed \
		> $@.part
	mv $@.part $@

%.coverage.exome.tsv: /data/christian/hdall/data/current/bam/%.merged.duplicate_marked.realigned.recalibrated.bam /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr \
		> $@.part
	mv $@.part $@

../reseq/cnv/%.coverage.reseq.panel.tsv: ~/hdall/data/reseq/bam/%.duplicate_marked.realigned.recalibrated.bam ~/hdall/data/reseq/enriched_regions.bed
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b ~/hdall/data/reseq/enriched_regions.bed \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed ~/hdall/data/reseq/enriched_regions.bed \
		> $@.part
	mv $@.part $@

%.coverage.mito.tsv: /data/christian/hdall/data/current/bam/%.merged.duplicate_marked.realigned.recalibrated.bam
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-r chrM \
		$< \
	> $@.part
	mv $@.part $@

cov-plot.%.exome.pdf: %_dia.coverage.exome.tsv %_rel.coverage.exome.tsv %_rem.coverage.exome.tsv ~/hdall/scripts/cnv/cov-plot-exome.R
	Rscript ~/hdall/scripts/cnv/cov-plot-exome.R --patient $* --diagnosis $(word 1,$^) --relapse $(word 2,$^) --remission $(word 3,$^) 2>&1
	mv cov-plot.$*.pdf.part $@

cov-plot.%.panel.pdf: %_dia.coverage.panel.tsv %_rel.coverage.panel.tsv %_rem.coverage.panel.tsv ~/hdall/scripts/cnv/cov-plot-panel.R
	Rscript ~/hdall/scripts/cnv/cov-plot-panel.R --patient $* --diagnosis $(word 1,$^) --relapse $(word 2,$^) --remission $(word 3,$^) 2>&1
	mv cov-plot.$*.reseq.pdf.part $@

../reseq/cnv/cov-plot.%.matched.reseq.panel.pdf: ../reseq/cnv/%_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/%_Relapse.coverage.reseq.panel.tsv ../reseq/cnv/%_Remission.coverage.reseq.panel.tsv ~/hdall/scripts/cnv/cov-plot-panel.R
	Rscript ~/hdall/scripts/cnv/cov-plot-panel.R --patient $* --diagnosis $(word 1,$^) --relapse $(word 2,$^) --remission $(word 3,$^) 2>&1
	mv cov-plot.$*.reseq.pdf.part $@
	
../reseq/cnv/cov-plot.%.dia.reseq.panel.pdf: ../reseq/cnv/%_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/%_Remission.coverage.reseq.panel.tsv ~/hdall/scripts/cnv/cov-plot-panel.R
	Rscript ~/hdall/scripts/cnv/cov-plot-panel.R --patient $* --diagnosis $(word 1,$^) --remission $(word 2,$^) 2>&1
	mv cov-plot.$*.reseq.pdf.part $@

../reseq/cnv/cov-plot.%.dia.nonrel.reseq.panel.pdf: ../reseq/cnv/%_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/1009302_Remission.coverage.reseq.panel.tsv ~/hdall/scripts/cnv/cov-plot-panel.R
	Rscript ~/hdall/scripts/cnv/cov-plot-panel.R --patient $* --diagnosis $(word 1,$^) --remission $(word 2,$^) 2>&1
	mv cov-plot.$*.reseq.pdf.part $@

../reseq/cnv/cov-plot.%.rel.reseq.panel.pdf: ../reseq/cnv/%_Relapse.coverage.reseq.panel.tsv ../reseq/cnv/%_Remission.coverage.reseq.panel.tsv ~/hdall/scripts/cnv/cov-plot-panel.R
	Rscript ~/hdall/scripts/cnv/cov-plot-panel.R --patient $* --relapse $(word 1,$^) --remission $(word 2,$^) 2>&1
	mv cov-plot.$*.reseq.pdf.part $@
	
../reseq/cnv/segmented_coverage/%.segmented-coverage.tsv: ~/hdall/data/reseq/bam/%.duplicate_marked.realigned.recalibrated.bam ~/hdall/scripts/cnv/get-segment-coverage.pl
	~/tools/samtools-0.1.19/samtools depth -Q 1 $< \
		| perl ~/hdall/scripts/cnv/get-segment-coverage.pl --sample $* --bin-size 250000 \
		> $@.part
	mv $@.part $@

../reseq/cnv/segmented_coverage/%_Diagnosis.segmented-coverage.pdf: ../reseq/cnv/segmented_coverage/%_Diagnosis.segmented-coverage.tsv ../reseq/cnv/segmented_coverage/%_Remission.segmented-coverage.tsv ~/hdall/scripts/cnv/plot-segment-coverage.R
	Rscript ~/hdall/scripts/cnv/plot-segment-coverage.R \
		--patient $*_dia \
		--tumor $(word 1,$^) \
		--normal $(word 2,$^) \
		--circos ~/hdall/results/cnv/circos/$*_dia.cnv.circos.tsv \
		--output $@.part
	mv $@.part $@

../reseq/cnv/segmented_coverage/%_Relapse.segmented-coverage.pdf: ../reseq/cnv/segmented_coverage/%_Relapse.segmented-coverage.tsv ../reseq/cnv/segmented_coverage/%_Remission.segmented-coverage.tsv ~/hdall/scripts/cnv/plot-segment-coverage.R
	Rscript ~/hdall/scripts/cnv/plot-segment-coverage.R \
		--patient $*_rel \
		--tumor $(word 1,$^) \
		--normal $(word 2,$^) \
		--circos ~/hdall/results/cnv/circos/$*_rel.cnv.circos.tsv \
		--output $@.part
	mv $@.part $@
	
../reseq/cnv/segmented_coverage/%_Diagnosis.nonrel.segmented-coverage.pdf: ../reseq/cnv/segmented_coverage/%_Diagnosis.segmented-coverage.tsv ../reseq/cnv/segmented_coverage/1009302_Remission.segmented-coverage.tsv ~/hdall/scripts/cnv/plot-segment-coverage.R
	Rscript ~/hdall/scripts/cnv/plot-segment-coverage.R \
		--patient $*_dia \
		--tumor $(word 1,$^) \
		--normal $(word 2,$^) \
		--circos ~/hdall/results/cnv/circos/$*_dia.nonrel.cnv.circos.tsv \
		--output $@.part
	mv $@.part $@
	
relapsing.rel.circos.png: $(foreach P, $(PATIENTS_RESEQ_MATCHED) $(PATIENTS_RESEQ_RELONLY), ~/hdall/results/cnv/circos/$P_rel.cnv.circos.tsv)
	 ~/tools/circos-0.64/bin/circos --conf ~/hdall/scripts/cnv/circos.conf
	 mv relapsing.rel.circos.part.png $@ 

#-----------------------------------------------------
# regional coverage plots
#-----------------------------------------------------
.PHONY: region-coverage
region-coverage: ../region-coverage/allpatients.IKZF2-chr2-211852462-217849831.pdf \
				 ../region-coverage/allpatients.IKZF1-chr7-49578046-51601231.pdf

reseq-region-coverage: ../reseq/region-coverage/allpatients.IKZF2-chr2-211852462-217849831.pdf \
					   ../reseq/region-coverage/allpatients.IKZF1-chr7-49578046-51601231.pdf

../region-coverage/allpatients.%.pdf: $(foreach P, $(PATIENTS_EXOME), ../region-coverage/patient$P.%.pdf)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@
	rm $^

../reseq/region-coverage/allpatients.%.pdf: $(foreach P, $(PATIENTS_EXOME), ../reseq/region-coverage/patient$P.%.exome.pdf) $(foreach P, $(PATIENTS_RESEQ_MATCHED), ../reseq/region-coverage/patient$P.%.reseq.pdf)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@
	rm $^

.SECONDEXPANSION: 
../region-coverage/patient%.pdf: $$(word 1, $$(subst ., , %))_dia.coverage.exome.tsv $$(word 1, $$(subst ., , %))_rel.coverage.exome.tsv $$(word 1, $$(subst ., , %))_rem.coverage.exome.tsv
	Rscript ~/p2ry8-crlf2/scripts/cov-plot-region.R \
		--patient $(word 1, $(subst ., , $*)) \
		--diagnosis $(word 1,$^) \
		--relapse $(word 2,$^) \
		--remission $(word 3,$^) \
		--output $@.part \
		--region-name $(word 1, $(subst -, , $(word 2, $(subst ., , $*)))) \
		--display-chrom $(word 2, $(subst -, , $*)) \
		--display-start $(word 3, $(subst -, , $*)) \
		--display-end $(word 4, $(subst -, , $*)) \
		$(if $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),--display-genes $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),)
	mv $@.part $@

../reseq/region-coverage/patient%.exome.pdf: ../reseq/cnv/$$(word 1, $$(subst ., , %))_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/$$(word 1, $$(subst ., , %))_Relapse.coverage.reseq.panel.tsv $$(word 1, $$(subst ., , %))_rem.coverage.panel.tsv
	Rscript ~/p2ry8-crlf2/scripts/cov-plot-region.R \
		--patient $(word 1, $(subst ., , $*)) \
		--diagnosis $(word 1,$^) \
		--relapse $(word 2,$^) \
		--remission $(word 3,$^) \
		--output $@.part \
		--region-name $(word 1, $(subst -, , $(word 2, $(subst ., , $*)))) \
		--display-chrom $(word 2, $(subst -, , $*)) \
		--display-start $(word 3, $(subst -, , $*)) \
		--display-end $(word 4, $(subst -, , $*)) \
		$(if $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),--display-genes $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),)
	mv $@.part $@

../reseq/region-coverage/patient%.reseq.pdf: ../reseq/cnv/$$(word 1, $$(subst ., , %))_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/$$(word 1, $$(subst ., , %))_Relapse.coverage.reseq.panel.tsv ../reseq/cnv/$$(word 1, $$(subst ., , %))_Remission.coverage.reseq.panel.tsv
	Rscript ~/p2ry8-crlf2/scripts/cov-plot-region.R \
		--patient $(word 1, $(subst ., , $*)) \
		--diagnosis $(word 1,$^) \
		--relapse $(word 2,$^) \
		--remission $(word 3,$^) \
		--output $@.part \
		--region-name $(word 1, $(subst -, , $(word 2, $(subst ., , $*)))) \
		--display-chrom $(word 2, $(subst -, , $*)) \
		--display-start $(word 3, $(subst -, , $*)) \
		--display-end $(word 4, $(subst -, , $*)) \
		$(if $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),--display-genes $(GENES_$(word 1, $(subst -, , $(word 2, $(subst ., , $*))))),)
	mv $@.part $@	