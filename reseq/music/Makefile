export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'
PATIENTS=1009302 1019357 1020540 1020583 1021247 1021392 1021631 1022914 1023056 1023392 1024518 1024543 1025108 1025409 1187 314 399 430 446 460 545 592 715 786 792 818 842 A B C D FB11 G HS6 K MB2 X 243 933 944 KA14651 KD20493 KE12025 MJ16441 NH17331 PJ13414 RD17412 RS13466 ST13892 ZA16211 Y
PATIENTS_REL_ONLY=1017005 1021865 1023545 AD15 BL16 BM18 CA18 DM1 FE1 FS1 GD18 GD1 HJ15 HJA15 HL1 KA17 KJ17 KL16 LB17 LM18 MJ1 ML10 MV16 NS18 PC16 RT15 RT16 SJM16 SKR1 SL1 SLM1 ST14 WA1 ZE13
PATIENTS_DIA_ONLY=1004564 1010661 1010781 1019964 1020076 1021087 1023338 1023616 1024589 1026233 1026662 B100 EF7 FB14 G44 HD7

all: mutation-relation-dia mutation-relation-rel dia/proximity.reseq.tsv rel/proximity.reseq.tsv 

panel-genes.roi: ../../panel-genes.tsv ../../kamilla/final-list-design-studio.tsv ~/hdall/scripts/reseq/music/get-roi-file.pl
	perl ~/hdall/scripts/reseq/music/get-roi-file.pl | sort -t '	' -k 1,1 -k 2,2g > $@.part
	mv $@.part $@

panel-genes.roi.gz.tbi: panel-genes.roi
	bgzip -c $^ > panel-genes.roi.gz
	~/tools/tabix-0.2.6/tabix panel-genes.roi.gz -s 1 -b 2 -e 3
	
rel/maf/%.reseq.maf: ../filtered_variants/%.reseq.dbsnp.snpeff.dbNSFP.filtered.vcf ../filtered_variants/%.reseq.indel.dbsnp.snpeff.dbNSFP.filtered.vcf panel-genes.roi.gz.tbi ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv ~/hdall/scripts/vcf2maf.pl
	mkdir -p rel/maf
	cat $(word 1,$^) | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi panel-genes.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.10 \
		2>&1 1>$@.part | $(LOG)
	cat $(word 2,$^) | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi panel-genes.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.10 \
		2>&1 1>>$@.part | $(LOG)
	mv $@.part $@

rel/maf/rem_rel.reseq.maf: $(foreach P, $(PATIENTS), rel/maf/$P_rem_rel.reseq.maf) $(foreach P, $(PATIENTS_REL_ONLY), rel/maf/$P_rem_rel.reseq.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>$@.part | $(LOG)
	cat $^ >> $@.part
	mv $@.part $@

rel/maf/rem_rel.reseq.matched.maf: $(foreach P, $(PATIENTS), rel/maf/$P_rem_rel.reseq.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>$@.part | $(LOG)
	cat $^ >> $@.part
	mv $@.part $@

dia/maf/rem_dia.reseq.maf: $(foreach P, $(PATIENTS), dia/maf/$P_rem_dia.reseq.maf) $(foreach P, $(PATIENTS_DIA_ONLY), dia/maf/$P_rem_dia.reseq.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>$@.part | $(LOG)
	cat $^ >> $@.part
	mv $@.part $@

dia/maf/%.reseq.maf: ../filtered_variants/%.reseq.dbsnp.snpeff.dbNSFP.filtered.vcf ../filtered_variants/%.reseq.indel.dbsnp.snpeff.dbNSFP.filtered.vcf panel-genes.roi.gz.tbi ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv ~/hdall/scripts/vcf2maf.pl
	mkdir -p dia/maf
	cat $(word 1,$^) | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi panel-genes.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.10 \
		2>&1 1>$@.part | $(LOG)
	cat $(word 2,$^) | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi panel-genes.roi.gz \
		--mapping-entrez ~/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.10 \
		2>&1 1>>$@.part | $(LOG)
	mv $@.part $@

#---------------------------
# mutation-relation relapse
#---------------------------
mutation-relation-rel: rel/mutation-relation.rel.af10.pathway.fine.reseq.tsv rel/mutation-relation.rel.af10.pathway.coarse.reseq.tsv rel/mutation-relation.recrel2.af10.reseq.tsv rel/mutation-relation.recrel2.combined-KRAS-PTPN11-NRAS.reseq.tsv rel/mutation-relation.recrel2.combined-KRAS-PTPN11.reseq.tsv rel/mutation-relation.recrel2.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.tsv rel/mutation-relation.recrel2.af10.reseq.matched.tsv

rel/wig-list.tsv:
	rm -f $@.part
	for P in $(PATIENTS) $(PATIENTS_REL_ONLY) ; do \
		echo "$$P""_rel	~/hdall/data/reseq/wig/$$P""_rel.merged.duplicate_marked.realigned.recalibrated.music.wig" >> $@.part; \
	done
	mv $@.part $@

rel/wig-list.matched.tsv:
	rm -f $@.part
	for P in $(PATIENTS) ; do \
		echo "$$P""_rel	~/hdall/data/reseq/wig/$$P""_rel.merged.duplicate_marked.realigned.recalibrated.music.wig" >> $@.part; \
	done
	mv $@.part $@

rel/maf/rem_rel.pathway.fine.reseq.maf: rel/maf/rem_rel.reseq.maf ~/hdall/scripts/reseq/aggregate-genes-to-pathways.pl
	cat $(word 1,$^) | perl ~/hdall/scripts/reseq/aggregate-genes-to-pathways.pl --level fine > $@.part
	mv $@.part $@
rel/mutation-relation.rel.af10.pathway.fine.reseq.gene-list.txt: rel/maf/rem_rel.pathway.fine.reseq.maf
	grep -vP "^(#|Hugo_Symbol)" $^ | cut -f 1 | sort | uniq > $@.part
	mv $@.part $@
rel/mutation-relation.rel.af10.pathway.fine.reseq.tsv: rel/mutation-relation.rel.af10.pathway.fine.reseq.gene-list.txt rel/wig-list.tsv rel/maf/rem_rel.pathway.fine.reseq.maf
	genome music mutation-relation \
		--bam-list=rel/wig-list.tsv \
		--maf-file=rel/maf/rem_rel.pathway.fine.reseq.maf \
		--permutations 1000 \
		--gene-list rel/mutation-relation.rel.af10.pathway.fine.reseq.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@

rel/maf/rem_rel.pathway.coarse.reseq.maf: rel/maf/rem_rel.reseq.maf ~/hdall/scripts/reseq/aggregate-genes-to-pathways.pl
	cat $(word 1,$^) | perl ~/hdall/scripts/reseq/aggregate-genes-to-pathways.pl --level coarse > $@.part
	mv $@.part $@
rel/mutation-relation.rel.af10.pathway.coarse.reseq.gene-list.txt: rel/maf/rem_rel.pathway.coarse.reseq.maf
	grep -vP "^(#|Hugo_Symbol)" $^ | cut -f 1 | sort | uniq > $@.part
	mv $@.part $@
rel/mutation-relation.rel.af10.pathway.coarse.reseq.tsv: rel/mutation-relation.rel.af10.pathway.coarse.reseq.gene-list.txt rel/wig-list.tsv rel/maf/rem_rel.pathway.coarse.reseq.maf
	genome music mutation-relation \
		--bam-list=rel/wig-list.tsv \
		--maf-file=rel/maf/rem_rel.pathway.coarse.reseq.maf \
		--permutations 1000 \
		--gene-list rel/mutation-relation.rel.af10.pathway.coarse.reseq.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@

rel/mutation-relation.recrel2.af10.gene-list.txt: ../gene-patient-matrix.reseq.tsv
	python ~/hdall/scripts/query-text.py \
		--file ../gene-patient-matrix.reseq.tsv \
		--sql "select gene,data.'max-af-rel-ns' AS m,data.'freq-rel-ns' AS f from data where m >= 10 and m <> '' and f >= 2" \
		> $@.part
	mv $@.part $@
rel/mutation-relation.recrel2.af10.reseq.tsv: rel/mutation-relation.recrel2.af10.gene-list.txt rel/wig-list.tsv rel/maf/rem_rel.reseq.maf
	genome music mutation-relation \
		--bam-list=rel/wig-list.tsv \
		--maf-file=rel/maf/rem_rel.reseq.maf \
		--permutations 1000 \
		--gene-list rel/mutation-relation.recrel2.af10.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@
rel/mutation-relation.recrel2.af10.reseq.matched.tsv: rel/mutation-relation.recrel2.af10.gene-list.txt rel/wig-list.matched.tsv rel/maf/rem_rel.reseq.matched.maf
	genome music mutation-relation \
		--bam-list=rel/wig-list.matched.tsv \
		--maf-file=rel/maf/rem_rel.reseq.matched.maf \
		--permutations 1000 \
		--gene-list rel/mutation-relation.recrel2.af10.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@


rel/maf/rem_rel.combined-KRAS-PTPN11-NRAS.reseq.maf: rel/maf/rem_rel.reseq.maf
	cat $^ | perl -ne 's/^(KRAS|PTPN11|NRAS)/KRAS-PTPN11-NRAS/; print $$_;' > $@.part
	mv $@.part $@
rel/mutation-relation.recrel2.af10.combined-KRAS-PTPN11-NRAS.reseq.gene-list.txt: rel/mutation-relation.recrel2.af10.gene-list.txt
	cat $^ | perl -ne 's/^(KRAS|PTPN11|NRAS)/KRAS-PTPN11-NRAS/; print $$_;' | cut -f 1 | sort | uniq > $@.part
	mv $@.part $@
rel/mutation-relation.recrel2.combined-KRAS-PTPN11-NRAS.reseq.tsv: rel/mutation-relation.recrel2.af10.combined-KRAS-PTPN11-NRAS.reseq.gene-list.txt rel/wig-list.tsv rel/maf/rem_rel.combined-KRAS-PTPN11-NRAS.reseq.maf
	genome music mutation-relation \
		--bam-list=rel/wig-list.tsv \
		--maf-file=rel/maf/rem_rel.combined-KRAS-PTPN11-NRAS.reseq.maf \
		--permutations 1000 \
		--gene-list rel/mutation-relation.recrel2.af10.combined-KRAS-PTPN11-NRAS.reseq.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@


rel/maf/rem_rel.combined-KRAS-PTPN11.reseq.maf: rel/maf/rem_rel.reseq.maf
	cat $^ | perl -ne 's/^(KRAS|PTPN11)/KRAS-PTPN11/; print $$_;' > $@.part
	mv $@.part $@
rel/mutation-relation.recrel2.af10.combined-KRAS-PTPN11.reseq.gene-list.txt: rel/mutation-relation.recrel2.af10.gene-list.txt
	cat $^ | perl -ne 's/^(KRAS|PTPN11)/KRAS-PTPN11/; print $$_;' | cut -f 1 | sort | uniq > $@.part
	mv $@.part $@
rel/mutation-relation.recrel2.combined-KRAS-PTPN11.reseq.tsv: rel/mutation-relation.recrel2.af10.combined-KRAS-PTPN11.reseq.gene-list.txt rel/wig-list.tsv rel/maf/rem_rel.combined-KRAS-PTPN11.reseq.maf
	genome music mutation-relation \
		--bam-list=rel/wig-list.tsv \
		--maf-file=rel/maf/rem_rel.combined-KRAS-PTPN11.reseq.maf \
		--permutations 1000 \
		--gene-list rel/mutation-relation.recrel2.af10.combined-KRAS-PTPN11.reseq.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@


rel/maf/rem_rel.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.maf: rel/maf/rem_rel.reseq.maf
	cat $^ | perl -ne 's/^(KRAS|PTPN11)/KRAS-PTPN11/; s/^(CREBBP|TRRAP)/CREBBP-TRRAP/; print $$_;' > $@.part
	mv $@.part $@
rel/mutation-relation.recrel2.af10.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.gene-list.txt: rel/mutation-relation.recrel2.af10.gene-list.txt
	cat $^ | perl -ne 's/^(KRAS|PTPN11)/KRAS-PTPN11/; s/^(CREBBP|TRRAP)/CREBBP-TRRAP/; print $$_;' | cut -f 1 | sort | uniq > $@.part
	mv $@.part $@
rel/mutation-relation.recrel2.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.tsv: rel/mutation-relation.recrel2.af10.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.gene-list.txt rel/wig-list.tsv rel/maf/rem_rel.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.maf
	genome music mutation-relation \
		--bam-list=rel/wig-list.tsv \
		--maf-file=rel/maf/rem_rel.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.maf \
		--permutations 1000 \
		--gene-list rel/mutation-relation.recrel2.af10.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@

#---------------------------
# mutation-relation diagnosis
#---------------------------
mutation-relation-dia: dia/mutation-relation.dia.af10.pathway.fine.reseq.tsv dia/mutation-relation.dia.af10.pathway.coarse.reseq.tsv dia/mutation-relation.recdia2.af10.reseq.tsv dia/mutation-relation.recdia2.combined-KRAS-PTPN11-NRAS.reseq.tsv dia/mutation-relation.recdia2.combined-KRAS-PTPN11.reseq.tsv dia/mutation-relation.recdia2.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.tsv

dia/maf/rem_dia.pathway.fine.reseq.maf: dia/maf/rem_dia.reseq.maf ~/hdall/scripts/reseq/aggregate-genes-to-pathways.pl
	cat $(word 1,$^) | perl ~/hdall/scripts/reseq/aggregate-genes-to-pathways.pl --level fine > $@.part
	mv $@.part $@
dia/mutation-relation.dia.af10.pathway.fine.reseq.gene-list.txt: dia/maf/rem_dia.pathway.fine.reseq.maf
	grep -vP "^(#|Hugo_Symbol)" $^ | cut -f 1 | sort | uniq > $@.part
	mv $@.part $@
dia/mutation-relation.dia.af10.pathway.fine.reseq.tsv: dia/mutation-relation.dia.af10.pathway.fine.reseq.gene-list.txt dia/wig-list.tsv dia/maf/rem_dia.pathway.fine.reseq.maf
	genome music mutation-relation \
		--bam-list=dia/wig-list.tsv \
		--maf-file=dia/maf/rem_dia.pathway.fine.reseq.maf \
		--permutations 1000 \
		--gene-list dia/mutation-relation.dia.af10.pathway.fine.reseq.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@

dia/maf/rem_dia.pathway.coarse.reseq.maf: dia/maf/rem_dia.reseq.maf ~/hdall/scripts/reseq/aggregate-genes-to-pathways.pl
	cat $(word 1,$^) | perl ~/hdall/scripts/reseq/aggregate-genes-to-pathways.pl --level coarse > $@.part
	mv $@.part $@
dia/mutation-relation.dia.af10.pathway.coarse.reseq.gene-list.txt: dia/maf/rem_dia.pathway.coarse.reseq.maf
	grep -vP "^(#|Hugo_Symbol)" $^ | cut -f 1 | sort | uniq > $@.part
	mv $@.part $@
dia/mutation-relation.dia.af10.pathway.coarse.reseq.tsv: dia/mutation-relation.dia.af10.pathway.coarse.reseq.gene-list.txt dia/wig-list.tsv dia/maf/rem_dia.pathway.coarse.reseq.maf
	genome music mutation-relation \
		--bam-list=dia/wig-list.tsv \
		--maf-file=dia/maf/rem_dia.pathway.coarse.reseq.maf \
		--permutations 1000 \
		--gene-list dia/mutation-relation.dia.af10.pathway.coarse.reseq.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@

dia/mutation-relation.recdia2.af10.gene-list.txt: ../gene-patient-matrix.reseq.tsv
	python ~/hdall/scripts/query-text.py \
		--file ../gene-patient-matrix.reseq.tsv \
		--sql "select gene,data.'max-af-dia-ns' AS m,data.'freq-dia-ns' AS f from data where m >= 10 and m <> '' and f >= 2" \
		> $@.part
	mv $@.part $@
dia/mutation-relation.recdia2.af10.reseq.tsv: dia/mutation-relation.recdia2.af10.gene-list.txt dia/wig-list.tsv dia/maf/rem_dia.reseq.maf
	genome music mutation-relation \
		--bam-list=dia/wig-list.tsv \
		--maf-file=dia/maf/rem_dia.reseq.maf \
		--permutations 1000 \
		--gene-list dia/mutation-relation.recdia2.af10.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@


dia/maf/rem_dia.combined-KRAS-PTPN11-NRAS.reseq.maf: dia/maf/rem_dia.reseq.maf
	cat $^ | perl -ne 's/^(KRAS|PTPN11|NRAS)/KRAS-PTPN11-NRAS/; print $$_;' > $@.part
	mv $@.part $@
dia/mutation-relation.recdia2.af10.combined-KRAS-PTPN11-NRAS.reseq.gene-list.txt: dia/mutation-relation.recdia2.af10.gene-list.txt
	cat $^ | perl -ne 's/^(KRAS|PTPN11|NRAS)/KRAS-PTPN11-NRAS/; print $$_;' | cut -f 1 | sort | uniq > $@.part
	mv $@.part $@
dia/mutation-relation.recdia2.combined-KRAS-PTPN11-NRAS.reseq.tsv: dia/mutation-relation.recdia2.af10.combined-KRAS-PTPN11-NRAS.reseq.gene-list.txt dia/wig-list.tsv dia/maf/rem_dia.combined-KRAS-PTPN11-NRAS.reseq.maf
	genome music mutation-relation \
		--bam-list=dia/wig-list.tsv \
		--maf-file=dia/maf/rem_dia.combined-KRAS-PTPN11-NRAS.reseq.maf \
		--permutations 1000 \
		--gene-list dia/mutation-relation.recdia2.af10.combined-KRAS-PTPN11-NRAS.reseq.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@


dia/maf/rem_dia.combined-KRAS-PTPN11.reseq.maf: dia/maf/rem_dia.reseq.maf
	cat $^ | perl -ne 's/^(KRAS|PTPN11)/KRAS-PTPN11/; print $$_;' > $@.part
	mv $@.part $@
dia/mutation-relation.recdia2.af10.combined-KRAS-PTPN11.reseq.gene-list.txt: dia/mutation-relation.recdia2.af10.gene-list.txt
	cat $^ | perl -ne 's/^(KRAS|PTPN11)/KRAS-PTPN11/; print $$_;' | cut -f 1 | sort | uniq > $@.part
	mv $@.part $@
dia/mutation-relation.recdia2.combined-KRAS-PTPN11.reseq.tsv: dia/mutation-relation.recdia2.af10.combined-KRAS-PTPN11.reseq.gene-list.txt dia/wig-list.tsv dia/maf/rem_dia.combined-KRAS-PTPN11.reseq.maf
	genome music mutation-relation \
		--bam-list=dia/wig-list.tsv \
		--maf-file=dia/maf/rem_dia.combined-KRAS-PTPN11.reseq.maf \
		--permutations 1000 \
		--gene-list dia/mutation-relation.recdia2.af10.combined-KRAS-PTPN11.reseq.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@

dia/maf/rem_dia.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.maf: dia/maf/rem_dia.reseq.maf
	cat $^ | perl -ne 's/^(KRAS|PTPN11)/KRAS-PTPN11/; s/^(CREBBP|TRRAP)/CREBBP-TRRAP/; print $$_;' > $@.part
	mv $@.part $@
dia/mutation-relation.recdia2.af10.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.gene-list.txt: dia/mutation-relation.recdia2.af10.gene-list.txt
	cat $^ | perl -ne 's/^(KRAS|PTPN11)/KRAS-PTPN11/; s/^(CREBBP|TRRAP)/CREBBP-TRRAP/; print $$_;' | cut -f 1 | sort | uniq > $@.part
	mv $@.part $@
dia/mutation-relation.recdia2.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.tsv: dia/mutation-relation.recdia2.af10.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.gene-list.txt dia/wig-list.tsv dia/maf/rem_dia.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.maf
	genome music mutation-relation \
		--bam-list=dia/wig-list.tsv \
		--maf-file=dia/maf/rem_dia.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.maf \
		--permutations 1000 \
		--gene-list dia/mutation-relation.recdia2.af10.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.gene-list.txt \
		--output-file $@.part \
		2>&1 | $(LOG)
	mv $@.part $@
	
#-----------------------------------------------------------------------
#--- GENOME MUSIC PROXIMITY
#-----------------------------------------------------------------------

dia/proximity.reseq.tsv: dia/maf/rem_dia.reseq.maf
	genome music proximity \
		--maf-file dia/maf/rem_dia.reseq.maf \
		--max-proximity 15 \
		--output-dir dia/ \
		2>&1 | tee -a make.log
	mv dia/proximity_report $@

rel/proximity.reseq.tsv: rel/maf/rem_rel.reseq.maf
	genome music proximity \
		--maf-file rel/maf/rem_rel.reseq.maf \
		--max-proximity 15 \
		--output-dir rel/ \
		2>&1 | tee -a make.log	
	mv rel/proximity_report $@

		