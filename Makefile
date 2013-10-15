export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash

all: filtered-variants.tsv filtered-variants.cosmic.tsv gene-patient-matrix.tsv gene-patient-matrix.high-af.tsv gene-patient-matrix.tier1.tsv cnv/gene-patient-matrix.cnv.tsv filtered-variants.cosmic.normaf.tsv lolliplot/lolliplot_CREBBP_NM_004380_both.svg ipa/mutated_relapse.tsv stats

stats: stats/mutations-per-patient-dia-vs-rel.pdf

nextera:
	cat kamilla/candidate\ genes\ for\ targeted\ sequencing.tsv | perl ~/hdall/scripts/get-nextera-exons.pl --density Standard > kamilla/nextera-exons.standard.csv
	cat kamilla/candidate\ genes\ for\ targeted\ sequencing.tsv | perl ~/hdall/scripts/get-nextera-exons.pl --density Dense > kamilla/nextera-exons.dense.csv
	cat kamilla/candidate\ genes\ for\ targeted\ sequencing.tsv | perl ~/hdall/scripts/get-nextera-exons.pl --density Standard --cds > kamilla/nextera-exons.standard.cds.csv
	cat kamilla/candidate\ genes\ for\ targeted\ sequencing.tsv | perl ~/hdall/scripts/get-nextera-exons.pl --density Dense --cds > kamilla/nextera-exons.dense.cds.csv

PATIENTS = 314 1021247 399 430 446 460 545 592 715 786 792 818 842 A B C D E X Y
filtered-variants.tsv:	$(foreach P, $(PATIENTS), filtered_variants/$P_rem_dia.snp.filtered.tsv) \
						$(foreach P, $(PATIENTS), filtered_variants/$P_rem_dia.indel.filtered.tsv) \
						$(foreach P, $(PATIENTS), filtered_variants/$P_rem_rel.snp.filtered.tsv) \
						$(foreach P, $(PATIENTS), filtered_variants/$P_rem_rel.indel.filtered.tsv) \
						~/hdall/scripts/filter-variants.pl 
	perl  ~/hdall/scripts/filter-variants.pl --header 2>&1 1>filtered-variants.tsv.part | tee -a make.log
	cat filtered_variants/*.filtered.tsv >> filtered-variants.tsv.part
	mv filtered-variants.tsv.part filtered-variants.tsv
	
filtered_variants/%.snp.filtered.tsv: ~/hdall/data/mutect_vcf/%_calls_snpeff_snpsift.vcf curated-recected-variants.tsv ~/hdall/scripts/filter-variants.pl
	perl ~/hdall/scripts/filter-variants.pl $* $< snp --vcf-out filtered_variants/$*.snp.filtered.vcf --rejected-variants-file curated-recected-variants.tsv \
		2>&1 1>$@.part | grep -v -P '(Leading or trailing space|variant.Format|on which the variant locates)' | tee -a make.log
	mv $@.part $@

filtered_variants/%.indel.filtered.tsv: ~/hdall/data/somatic_indel_vcf/%_snpeff.vcf curated-recected-variants.tsv ~/hdall/scripts/filter-variants.pl
	perl ~/hdall/scripts/filter-variants.pl $* $< indel --vcf-out filtered_variants/$*.indel.filtered.vcf --rejected-variants-file curated-recected-variants.tsv \
		2>&1 1>$@.part | grep -v -P '(Leading or trailing space|variant.Format|on which the variant locates)' | tee -a make.log
	mv $@.part $@	

filtered-variants.cosmic.tsv: filtered-variants.tsv ~/hdall/data/cosmic/v65/CosmicMutantExport_v65_280513.tsv ~/hdall/scripts/annotate-cosmic.pl
	cat ~/hdall/results/filtered-variants.tsv | perl ~/hdall/scripts/annotate-cosmic.pl \
		--cosmic-mutation-file ~/hdall/data/cosmic/v65/CosmicMutantExport_v65_280513.tsv \
		--only-confirmed \
		2>&1 1>$@.part | tee -a make.log
	mv $@.part $@ 

filtered-variants.cosmic.normaf.tsv: filtered-variants.cosmic.tsv cnv/hdall.cnv.tsv ~/hdall/scripts/normalize-af.pl
	cat filtered-variants.cosmic.tsv | perl ~/hdall/scripts/normalize-af.pl \
		--cnv-file cnv/hdall.cnv.tsv \
		2>&1 1>$@.part | tee -a make.log
	mv $@.part $@ 

impacted-genes-list.tsv: filtered-variants.tsv ~/hdall/scripts/impacted-genes.pl
	cat filtered-variants.tsv | perl ~/hdall/scripts/impacted-genes.pl \
		2>&1 1>impacted-genes-list.tsv.part | tee -a make.log
	mv impacted-genes-list.tsv.part impacted-genes-list.tsv

impacted-genes-list.high-af.tsv: filtered-variants.tsv ~/hdall/scripts/impacted-genes.pl
	cat filtered-variants.tsv \
		| perl -ne 'print $$_ if ((split/\t/)[19] >= 0.20)' \
		| perl ~/hdall/scripts/impacted-genes.pl \
		2>&1 1>impacted-genes-list.high-af.tsv.part | tee -a make.log
	mv impacted-genes-list.high-af.tsv.part impacted-genes-list.high-af.tsv

gene-patient-matrix.tsv: impacted-genes-list.tsv ~/hdall/scripts/get-gene-patient-matrix.pl
	cat impacted-genes-list.tsv | perl ~/hdall/scripts/get-gene-patient-matrix.pl --mut-details \
		2>&1 1>gene-patient-matrix.tsv.part | tee -a make.log
	mv gene-patient-matrix.tsv.part gene-patient-matrix.tsv

gene-patient-matrix.high-af.tsv: impacted-genes-list.high-af.tsv ~/hdall/scripts/get-gene-patient-matrix.pl
	cat impacted-genes-list.high-af.tsv | perl ~/hdall/scripts/get-gene-patient-matrix.pl --mut-details \
		2>&1 1>gene-patient-matrix.high-af.tsv.part | tee -a make.log
	mv gene-patient-matrix.high-af.tsv.part gene-patient-matrix.high-af.tsv

impacted-genes-list.tier1.tsv: filtered-variants.tsv ~/hdall/scripts/impacted-genes.pl
	cat filtered-variants.tsv \
		| perl -ne 'print $$_ if (/(patient|HIGH|MODERATE)/ and (split/\t/)[19] >= 0.20)' \
		| perl ~/hdall/scripts/impacted-genes.pl \
		2>&1 1>impacted-genes-list.tier1.tsv.part | tee -a make.log
	mv impacted-genes-list.tier1.tsv.part impacted-genes-list.tier1.tsv

gene-patient-matrix.tier1.tsv: impacted-genes-list.tier1.tsv ~/hdall/scripts/get-gene-patient-matrix.pl
	cat impacted-genes-list.tier1.tsv | perl ~/hdall/scripts/get-gene-patient-matrix.pl --mut-details \
		2>&1 1>gene-patient-matrix.tier1.tsv.part | tee -a make.log
	mv gene-patient-matrix.tier1.tsv.part gene-patient-matrix.tier1.tsv

cnv/hdall.cnv.tsv: cnv/CNAsallpatients.tsv cnv/hyperdiploid_CytoHDarray.tsv ~/hdall/scripts/cnv/parse-cnv-data.pl
	cat cnv/CNAsallpatients.tsv | perl ~/hdall/scripts/cnv/parse-cnv-data.pl --format andrea --header \
		2>&1 1>cnv/hdall.cnv.tsv.part | tee -a make.log 
	cat cnv/hyperdiploid_CytoHDarray.tsv | perl ~/hdall/scripts/cnv/parse-cnv-data.pl --format maria \
		2>&1 1>>cnv/hdall.cnv.tsv.part | tee -a make.log 
	mv cnv/hdall.cnv.tsv.part cnv/hdall.cnv.tsv

cnv/impacted-genes-list.cnv.tsv: cnv/hdall.cnv.tsv ~/hdall/scripts/cnv/impacted-genes-cnv.pl
	cat cnv/hdall.cnv.tsv | perl ~/hdall/scripts/cnv/impacted-genes-cnv.pl \
		--max-genes 99999 \
		2>&1 1>cnv/impacted-genes-list.cnv.tsv.part | tee -a make.log 
	mv cnv/impacted-genes-list.cnv.tsv.part cnv/impacted-genes-list.cnv.tsv

cnv/gene-patient-matrix.cnv.tsv: cnv/impacted-genes-list.cnv.tsv ~/hdall/scripts/get-gene-patient-matrix.cnv.pl
	cat cnv/impacted-genes-list.cnv.tsv | perl ~/hdall/scripts/get-gene-patient-matrix.cnv.pl \
		--max-genes 5 \
		2>&1 1>cnv/gene-patient-matrix.cnv.tsv.part | tee -a make.log
	mv cnv/gene-patient-matrix.cnv.tsv.part cnv/gene-patient-matrix.cnv.tsv

stats/variants-per-chrom-and-ploidy.pdf: filtered-variants.cosmic.normaf.tsv ~/hdall/scripts/somatic-variants-stat.R
	R --no-save --quiet --slave -f ~/hdall/scripts/somatic-variants-stat.R \
		2>&1 | tee -a make.log

clonal-analysis/allelic-freq-prob.distributions.pdf: ~/hdall/scripts/clonal-analysis/estimate-af-dist.R
	R --no-save --quiet --slave -f ~/hdall/scripts/clonal-analysis/estimate-af-dist.R \
		2>&1 | tee -a make.log

lolliplot/pfam-regions.filtered.tsv: ~/hdall/data/pfam-27.0/pfamA.txt ~/hdall/data/pfam-27.0/Pfam-A.regions.tsv.gz ~/hdall/scripts/lolliplot/filter-and-annotate-pfam-regions.pl
	perl ~/hdall/scripts/lolliplot/filter-and-annotate-pfam-regions.pl > lolliplot/pfam-regions.filtered.tsv

lolliplot/lolliplot_CREBBP_NM_004380_both.svg: id-mappings.tsv lolliplot/pfam-regions.filtered.tsv filtered-variants.cosmic.tsv ~/hdall/scripts/lolliplot/lolliplot.pl
	rm lolliplot/*.svg 
	perl ~/hdall/scripts/lolliplot/lolliplot.pl \
		2>&1 | tee -a make.log

ipa/mutated_relapse.tsv: music/rel-high-af/rem_rel.maf music/rel-high-af/smg.tsv
	cut -f 1,9 music/rel-high-af/rem_rel.maf | grep -P '(Frame_Shift|In_Frame|Missense|Nonsense|Splice_Site)' | cut -f 1 | sort | uniq | grep -f - music/rel-high-af/smg.tsv | cut -f 1,9 | sort -k 2g > ipa/mutated_relapse.tsv.tmp
	grep -vP "^(MESP2|TTN|TBP|CSMD3|DNAH5|RYR1|RYR2|RYR3|DNAH1|DNAH8|DNAH9|MUC2|MUC16|MUC12|MUC5B|OR5H2|OR5H2|OR11H4|OR6F1|OR52R1|OR2T12|OR6V1|OR51I2|OR5I1|OR9Q1|OR9Q1|PDZD7)\t" ipa/mutated_relapse.tsv.tmp > ipa/mutated_relapse.tsv.part
	rm ipa/mutated_relapse.tsv.tmp
	mv ipa/mutated_relapse.tsv.part ipa/mutated_relapse.tsv

ipa/mutated_relapse.noKRASpatients.tsv: music/rel-high-af/rem_rel.maf music/rel-high-af/smg.tsv
	cut -f 1,9,16 music/rel-high-af/rem_rel.maf | grep -P '(1021247|818|842|C|786|X|592|D|399|446|314|792)_rel' | grep -P '(Frame_Shift|In_Frame|Missense|Nonsense|Splice_Site)' | cut -f 1 | sort | uniq | grep -f - music/rel-high-af/smg.tsv | cut -f 1,9 | sort -k 2g > ipa/mutated_relapse.noKRASpatients.tsv.tmp
	grep -vP "^(MESP2|TTN|TBP|CSMD3|DNAH5|RYR1|RYR2|RYR3|DNAH1|DNAH8|DNAH9|MUC2|MUC16|MUC12|MUC5B|OR5H2|OR5H2|OR11H4|OR6F1|OR52R1|OR2T12|OR6V1|OR51I2|OR5I1|OR9Q1|OR9Q1|PDZD7)\t" ipa/mutated_relapse.noKRASpatients.tsv.tmp > ipa/mutated_relapse.noKRASpatients.tsv.part
	rm ipa/mutated_relapse.noKRASpatients.tsv.tmp
	mv ipa/mutated_relapse.noKRASpatients.tsv.part ipa/mutated_relapse.noKRASpatients.tsv

stats/mutations-per-patient-dia-vs-rel.pdf: filtered-variants.cosmic.normaf.tsv 
	R --no-save --quiet --slave -f ~/hdall/scripts/stats/mutations-per-patient.R 2>&1 | tee -a make.log
