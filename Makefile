export SHELLOPTS:=errexit:pipefail

all: gene-patient-matrix.tsv gene-patient-matrix.tier1.tsv

impacted-genes-list.tsv: filtered-variants.tsv ~/git/hdall/impacted-genes.pl
	cat filtered-variants.tsv | perl ~/git/hdall/impacted-genes.pl \
		2>&1 1>impacted-genes-list.tsv.part | tee -a make.log
	mv impacted-genes-list.tsv.part impacted-genes-list.tsv

gene-patient-matrix.tsv: impacted-genes-list.tsv ~/git/hdall/get-gene-patient-matrix.pl
	cat impacted-genes-list.tsv | perl ~/git/hdall/get-gene-patient-matrix.pl --mut-count \
		2>&1 1>gene-patient-matrix.tsv.part | tee -a make.log
	mv gene-patient-matrix.tsv.part gene-patient-matrix.tsv

	cat impacted-genes-list.tsv | perl ~/git/hdall/get-gene-patient-matrix.pl --mut-max-freq \
		2>&1 1>gene-patient-matrix.maxfreq.tsv.part | tee -a make.log
	mv gene-patient-matrix.maxfreq.tsv.part gene-patient-matrix.maxfreq.tsv
		
	cat impacted-genes-list.tsv | perl ~/git/hdall/get-gene-patient-matrix.pl --mut-details \
		2>&1 1>gene-patient-matrix.details.tsv.part | tee -a make.log
	mv gene-patient-matrix.details.tsv.part gene-patient-matrix.details.tsv

impacted-genes-list.tier1.tsv: filtered-variants.tsv ~/git/hdall/impacted-genes.pl
	cat filtered-variants.tsv \
		| perl -ne 'print $$_ if (/(patient|HIGH|MODERATE)/ and (split/\t/)[18] >= 0.25)' \
		| perl ~/git/hdall/impacted-genes.pl \
		2>&1 1>impacted-genes-list.tier1.tsv.part | tee -a make.log
	mv impacted-genes-list.tier1.tsv.part impacted-genes-list.tier1.tsv

gene-patient-matrix.tier1.tsv: impacted-genes-list.tier1.tsv ~/git/hdall/get-gene-patient-matrix.pl
	cat impacted-genes-list.tier1.tsv | perl ~/git/hdall/get-gene-patient-matrix.pl --mut-details \
		2>&1 1>gene-patient-matrix.tier1.tsv.part | tee -a make.log
	mv gene-patient-matrix.tier1.tsv.part gene-patient-matrix.tier1.tsv
