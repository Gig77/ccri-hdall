export SHELLOPTS:=errexit:pipefail

all: filtered-variants.tsv gene-patient-matrix.tsv gene-patient-matrix.tier1.tsv

PATIENTS = 314 1021247 399 430 446 460 545 592 715 786 792 818 842 A B C D E X Y
filtered-variants.tsv:	$(foreach P, $(PATIENTS), filtered_variants/$P_rem_dia.snp.filtered.tsv) \
						$(foreach P, $(PATIENTS), filtered_variants/$P_rem_dia.indel.filtered.tsv) \
						$(foreach P, $(PATIENTS), filtered_variants/$P_rem_rel.snp.filtered.tsv) \
						$(foreach P, $(PATIENTS), filtered_variants/$P_rem_rel.indel.filtered.tsv)
	perl  ~/hdall/scripts/filter-variants.pl --header 2>&1 1>filtered-variants.tsv.part | tee -a make.log
	cat filtered_variants/*.filtered.tsv >> filtered-variants.tsv.part
	mv filtered-variants.tsv.part filtered-variants.tsv
	
filtered_variants/%.snp.filtered.tsv: ~/hdall/data/mutect_vcf/%_calls_snpeff.vcf
	perl ~/hdall/scripts/filter-variants.pl $* $< snp --vcf-out filtered_variants/$*.snp.filtered.vcf \
		2>&1 1>$@ | grep -v -P '(Leading or trailing space|variant.Format)' | tee -a make.log

filtered_variants/%.indel.filtered.tsv: ~/hdall/data/somatic_indel_vcf/%_snpeff.vcf
	perl ~/hdall/scripts/filter-variants.pl $* $< indel --vcf-out filtered_variants/$*.indel.filtered.vcf \
		2>&1 1>$@ | grep -v -P '(Leading or trailing space|variant.Format)' | tee -a make.log

impacted-genes-list.tsv: filtered-variants.tsv ~/hdall/scripts/impacted-genes.pl
	cat filtered-variants.tsv | perl ~/hdall/scripts/impacted-genes.pl \
		2>&1 1>impacted-genes-list.tsv.part | tee -a make.log
	mv impacted-genes-list.tsv.part impacted-genes-list.tsv

gene-patient-matrix.tsv: impacted-genes-list.tsv ~/hdall/scripts/get-gene-patient-matrix.pl
	cat impacted-genes-list.tsv | perl ~/hdall/scripts/get-gene-patient-matrix.pl --mut-count \
		2>&1 1>gene-patient-matrix.tsv.part | tee -a make.log
	mv gene-patient-matrix.tsv.part gene-patient-matrix.tsv

	cat impacted-genes-list.tsv | perl ~/hdall/scripts/get-gene-patient-matrix.pl --mut-max-freq \
		2>&1 1>gene-patient-matrix.maxfreq.tsv.part | tee -a make.log
	mv gene-patient-matrix.maxfreq.tsv.part gene-patient-matrix.maxfreq.tsv
		
	cat impacted-genes-list.tsv | perl ~/hdall/scripts/get-gene-patient-matrix.pl --mut-details \
		2>&1 1>gene-patient-matrix.details.tsv.part | tee -a make.log
	mv gene-patient-matrix.details.tsv.part gene-patient-matrix.details.tsv

impacted-genes-list.tier1.tsv: filtered-variants.tsv ~/hdall/scripts/impacted-genes.pl
	cat filtered-variants.tsv \
		| perl -ne 'print $$_ if (/(patient|HIGH|MODERATE)/ and (split/\t/)[18] >= 0.25)' \
		| perl ~/hdall/scripts/impacted-genes.pl \
		2>&1 1>impacted-genes-list.tier1.tsv.part | tee -a make.log
	mv impacted-genes-list.tier1.tsv.part impacted-genes-list.tier1.tsv

gene-patient-matrix.tier1.tsv: impacted-genes-list.tier1.tsv ~/hdall/scripts/get-gene-patient-matrix.pl
	cat impacted-genes-list.tier1.tsv | perl ~/hdall/scripts/get-gene-patient-matrix.pl --mut-details \
		2>&1 1>gene-patient-matrix.tier1.tsv.part | tee -a make.log
	mv gene-patient-matrix.tier1.tsv.part gene-patient-matrix.tier1.tsv
