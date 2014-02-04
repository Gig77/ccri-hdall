export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

#PATIENTS_EXOME = 314 399 430 446 460 545 592 715 786 792 818 842 A B C D E X Y
#PATIENTS_RESEQ_MATCHED = 1009302 1019357 1020540 1020583 1021392 1021631 1022914 1023056 1023392 1024518 1024543 1025108 1025409 1187 FB11 G HS6 K MB2 243 933 944 KA14651 KD20493 KE12025 MJ16441 NH17331 PJ13414 RD17412 RS13466 ST13892 ZA16211
#PATIENTS_RESEQ_RELONLY = 1017005 1021865 1023545 AD15 BL16 BM18 CA18 DM1 FE1 FS1 GD18 GD1 HJ15 HJA15 HL1 KA17 KJ17 KL16 LB17 LM18 MJ1 ML10 MV16 NS18 PC16 RT15 RT16 SJM16 SKR1 SL1 SLM1 ST14 WA1 ZE13
#PATIENTS_RESEQ_DIAONLY = 1004564 1010661 1010781 1019964 1020076 1021087 1023338 1023616 1024589 1026233 1026662 B100 EF7 FB14 G44 HD7
PATIENTS_EXOME = 
PATIENTS_RESEQ_MATCHED = 
PATIENTS_RESEQ_RELONLY = 
PATIENTS_RESEQ_DIAONLY = 1010661 
NOT_WORKING = 1021247


all: $(foreach P, $(PATIENTS_EXOME), $P_rem.coverage.exome.tsv $P_dia.coverage.exome.tsv $P_rel.coverage.exome.tsv $P_rem.coverage.panel.tsv $P_dia.coverage.panel.tsv $P_rel.coverage.panel.tsv cov-plot.$P.exome.pdf cov-plot.$P.panel.pdf) \
	 $(foreach P, $(PATIENTS_RESEQ_MATCHED), ../reseq/cnv/$P_Remission.coverage.reseq.panel.tsv ../reseq/cnv/$P_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/$P_Relapse.coverage.reseq.panel.tsv ../reseq/cnv/cov-plot.$P.matched.reseq.panel.pdf) \
	 $(foreach P, $(PATIENTS_RESEQ_RELONLY), ../reseq/cnv/$P_Remission.coverage.reseq.panel.tsv ../reseq/cnv/$P_Relapse.coverage.reseq.panel.tsv ../reseq/cnv/cov-plot.$P.rel.reseq.panel.pdf) \
	 $(foreach P, $(PATIENTS_RESEQ_DIAONLY), ../reseq/cnv/$P_Remission.coverage.reseq.panel.tsv ../reseq/cnv/$P_Diagnosis.coverage.reseq.panel.tsv ../reseq/cnv/cov-plot.$P.dia.reseq.panel.pdf)

%.coverage.panel.tsv: /data/christian/hdall/data/current/bam/%.merged.duplicate_marked.realigned.recalibrated.bam /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr ~/git/hdall/cnv/get-exon-coverage.pl
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b ~/hdall/data/reseq/enriched_regions.bed \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed ~/hdall/data/reseq/enriched_regions.bed \
		> $@.part
	mv $@.part $@

%.coverage.exome.tsv: /data/christian/hdall/data/current/bam/%.merged.duplicate_marked.realigned.recalibrated.bam /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr ~/git/hdall/cnv/get-exon-coverage.pl
	~/tools/samtools-0.1.19/samtools depth \
		-Q 1 \
		-b /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr \
		$< \
	| perl ~/git/hdall/cnv/get-exon-coverage.pl \
		--exon-bed /data/christian/generic/data/current/illumina/truseq_exome_targeted_regions.hg19.bed.chr \
		> $@.part
	mv $@.part $@

../reseq/cnv/%.coverage.reseq.panel.tsv: ~/hdall/data/reseq/bam/%.duplicate_marked.realigned.recalibrated.bam ~/hdall/data/reseq/enriched_regions.bed ~/git/hdall/cnv/get-exon-coverage.pl
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

../reseq/cnv/cov-plot.%.rel.reseq.panel.pdf: ../reseq/cnv/%_Relapse.coverage.reseq.panel.tsv ../reseq/cnv/%_Remission.coverage.reseq.panel.tsv ~/hdall/scripts/cnv/cov-plot-panel.R
	Rscript ~/hdall/scripts/cnv/cov-plot-panel.R --patient $* --relapse $(word 1,$^) --remission $(word 2,$^) 2>&1
	mv cov-plot.$*.reseq.pdf.part $@
	