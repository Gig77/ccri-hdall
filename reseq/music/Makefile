export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'
PATIENTS=1009302 1019357 1020540 1020583 1021247 1021392 1021631 1022914 1023056 1023392 1024518 1024543 1025108 1025409 1187 314 399 430 446 460 545 592 715 786 792 818 842 A B C D FB11 G HS6 K MB2 X

all: mutation-relation-dia mutation-relation-rel 

panel-genes.roi: ../../panel-genes.tsv ../../kamilla/final-list-design-studio.tsv ~/hdall/scripts/reseq/music/get-roi-file.pl
	perl ~/hdall/scripts/reseq/music/get-roi-file.pl | sort -t '	' -k 1,1 -k 2,2g > $@.part
	mv $@.part $@

panel-genes.roi.gz.tbi: panel-genes.roi
	bgzip -c $^ > panel-genes.roi.gz
	~/tools/tabix-0.2.6/tabix panel-genes.roi.gz -s 1 -b 2 -e 3
	
rel/maf/%.reseq.maf: ../filtered_variants/%.reseq.dbsnp.snpeff.dbNSFP.filtered.vcf ../filtered_variants/%.reseq.indel.dbsnp.snpeff.dbNSFP.filtered.vcf panel-genes.roi.gz.tbi ~/hdall/results/gene-id-mapping.biomart-0.7.tsv ~/hdall/scripts/vcf2maf.pl
	mkdir -p rel/maf
	cat $(word 1,$^) | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi panel-genes.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.10 \
		2>&1 1>$@.part | $(LOG)
	cat $(word 2,$^) | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi panel-genes.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.10 \
		2>&1 1>>$@.part | $(LOG)
	mv $@.part $@

rel/maf/rem_rel.reseq.maf: $(foreach P, $(PATIENTS), rel/maf/$P_rem_rel.reseq.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>$@.part | $(LOG)
	cat rel/maf/*rem_rel.reseq.maf >> $@.part
	mv $@.part $@

dia/maf/rem_dia.reseq.maf: $(foreach P, $(PATIENTS), dia/maf/$P_rem_dia.reseq.maf)
	perl ~/hdall/scripts/vcf2maf.pl --header 2>&1 1>$@.part | $(LOG)
	cat dia/maf/*rem_dia.reseq.maf >> $@.part
	mv $@.part $@

dia/maf/%.reseq.maf: ../filtered_variants/%.reseq.dbsnp.snpeff.dbNSFP.filtered.vcf ../filtered_variants/%.reseq.indel.dbsnp.snpeff.dbNSFP.filtered.vcf panel-genes.roi.gz.tbi ~/hdall/results/gene-id-mapping.biomart-0.7.tsv ~/hdall/scripts/vcf2maf.pl
	mkdir -p dia/maf
	cat $(word 1,$^) | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi panel-genes.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.10 \
		2>&1 1>$@.part | $(LOG)
	cat $(word 2,$^) | perl ~/hdall/scripts/vcf2maf.pl \
		--sample $* \
		--music-roi panel-genes.roi.gz \
		--mapping-entrez ~/hdall/results/gene-id-mapping.biomart-0.7.tsv \
		--min-af 0.10 \
		2>&1 1>>$@.part | $(LOG)
	mv $@.part $@

#---------------------------
# mutation-relation relapse
#---------------------------
mutation-relation-rel: rel/mutation-relation.rel.af10.pathway.fine.reseq.tsv rel/mutation-relation.rel.af10.pathway.coarse.reseq.tsv rel/mutation-relation.recrel2.af10.reseq.tsv rel/mutation-relation.recrel2.combined-KRAS-PTPN11-NRAS.reseq.tsv rel/mutation-relation.recrel2.combined-KRAS-PTPN11.reseq.tsv rel/mutation-relation.recrel2.combined-CREBBP-TRRAP-KRAS-PTPN11.reseq.tsv

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
		