# download required tables from UCSC genome browser:
# mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select * from knownGene' > ~/generic/data/hg19/hg19.knownGene.txt
# mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select * from knownCanonical' > ~/generic/data/hg19/hg19.knownCanonical.txt
# mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'select * from kgXref' > ~/generic/data/hg19/hg19.kgXref.txt
# perl ~/hdall/scripts/get-id-mapping.pl > ~/hdall/results/id-mappings.tsv

SCRIPT=~/git/hdall/filter-variants.pl
VCF_DIR=/home/STANNANET/christian.frech/hdall/data/mutect_vcf
VCF_DIR_INDEL=/home/STANNANET/christian.frech/hdall/data/somatic_indel_vcf
VCF_OUT_FILTERED=/home/STANNANET/christian.frech/hdall/results/filtered_variants
VCF_OUT_MERGED=/home/STANNANET/christian.frech/hdall/results/merged_vcf
RESULT_FILE=/home/STANNANET/christian.frech/hdall/results/filtered-variants.tsv
LOG_FILE=/home/STANNANET/christian.frech/hdall/results/filtered-variants.log

perl $SCRIPT 314 rem_dia $VCF_DIR/314_rem_dia_calls_snpeff.vcf 314_rem 314_dia snp --header --vcf-out $VCF_OUT_FILTERED/314_rem_dia.snp.filtered.vcf >$RESULT_FILE 2>$LOG_FILE
perl $SCRIPT 314 rem_dia $VCF_DIR_INDEL/314_rem_dia_snpeff.vcf 314_rem 314_dia indel --vcf-out $VCF_OUT_FILTERED/314_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 314 rem_rel $VCF_DIR/314_rem_rel_calls_snpeff.vcf 314_rem 314_rel snp --vcf-out $VCF_OUT_FILTERED/314_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 314 rem_rel $VCF_DIR_INDEL/314_rem_rel_snpeff.vcf 314_rem 314_rel indel --vcf-out $VCF_OUT_FILTERED/314_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

#NOTE: merging/concatenating snp and indel VCFs is not that straightforward because they differ in their header, and vcftools does not deal with that properly;
#so we keep them separate for now
#~/tools/vcftools_0.1.10/bin/vcf-concat $VCF_OUT_FILTERED/314_rem_dia_calls_snpeff.snp.filtered.vcf $VCF_OUT_FILTERED/314_rem_dia_snpeff.indel.filtered.vcf 2>>$LOG_FILE | ~/tools/vcftools_0.1.10/bin/vcf-sort > $VCF_OUT_MERGED/314_rem_dia.filtered.merged.vcf   
#~/tools/vcftools_0.1.10/bin/vcf-concat $VCF_OUT_FILTERED/314_rem_rel_calls_snpeff.snp.filtered.vcf $VCF_OUT_FILTERED/314_rem_rel_snpeff.indel.filtered.vcf 2>>$LOG_FILE | ~/tools/vcftools_0.1.10/bin/vcf-sort > $VCF_OUT_MERGED/314_rem_rel.filtered.merged.vcf   

perl $SCRIPT 399 rem_dia $VCF_DIR/399_rem_dia_calls_snpeff.vcf 399_rem 399_dia snp --vcf-out $VCF_OUT_FILTERED/399_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 399 rem_dia $VCF_DIR_INDEL/399_rem_dia_snpeff.vcf 399_rem 399_dia indel --vcf-out $VCF_OUT_FILTERED/399_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 399 rem_rel $VCF_DIR/399_rem_rel_calls_snpeff.vcf 399_rem 399_rel snp --vcf-out $VCF_OUT_FILTERED/399_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 399 rem_rel $VCF_DIR_INDEL/399_rem_rel_snpeff.vcf 399_rem 399_rel indel --vcf-out $VCF_OUT_FILTERED/399_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 430 rem_dia $VCF_DIR/430_rem_dia_calls_snpeff.vcf 430_rem 430_dia snp --vcf-out $VCF_OUT_FILTERED/430_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 430 rem_dia $VCF_DIR_INDEL/430_rem_dia_snpeff.vcf 430_rem 430_dia indel --vcf-out $VCF_OUT_FILTERED/430_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 430 rem_rel $VCF_DIR/430_rem_rel_calls_snpeff.vcf 430_rem 430_rel snp --vcf-out $VCF_OUT_FILTERED/430_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 430 rem_rel $VCF_DIR_INDEL/430_rem_rel_snpeff.vcf 430_rem 430_rel indel --vcf-out $VCF_OUT_FILTERED/430_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 446 rem_dia $VCF_DIR/446_rem_dia_calls_snpeff.vcf 446_rem 446_dia snp --vcf-out $VCF_OUT_FILTERED/446_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 446 rem_dia $VCF_DIR_INDEL/446_rem_dia_snpeff.vcf 446_rem 446_dia indel --vcf-out $VCF_OUT_FILTERED/446_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 446 rem_rel $VCF_DIR/446_rem_rel_calls_snpeff.vcf 446_rem 446_rel snp --vcf-out $VCF_OUT_FILTERED/446_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 446 rem_rel $VCF_DIR_INDEL/446_rem_rel_snpeff.vcf 446_rem 446_rel indel --vcf-out $VCF_OUT_FILTERED/446_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 460 rem_dia $VCF_DIR/460_rem_dia_calls_snpeff.vcf 460_rem 460_dia snp --vcf-out $VCF_OUT_FILTERED/460_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 460 rem_dia $VCF_DIR_INDEL/460_rem_dia_snpeff.vcf 460_rem 460_dia indel --vcf-out $VCF_OUT_FILTERED/460_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 460 rem_rel $VCF_DIR/460_rem_rel_calls_snpeff.vcf 460_rem 460_rel snp --vcf-out $VCF_OUT_FILTERED/460_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 460 rem_rel $VCF_DIR_INDEL/460_rem_rel_snpeff.vcf 460_rem 460_rel indel --vcf-out $VCF_OUT_FILTERED/460_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 545 rem_dia $VCF_DIR/545_rem_dia_calls_snpeff.vcf 545_rem 545_dia snp --vcf-out $VCF_OUT_FILTERED/545_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 545 rem_dia $VCF_DIR_INDEL/545_rem_dia_snpeff.vcf 545_rem 545_dia indel --vcf-out $VCF_OUT_FILTERED/545_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 545 rem_rel $VCF_DIR/545_rem_rel_calls_snpeff.vcf 545_rem 545_rel snp --vcf-out $VCF_OUT_FILTERED/545_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 545 rem_rel $VCF_DIR_INDEL/545_rem_rel_snpeff.vcf 545_rem 545_rel indel --vcf-out $VCF_OUT_FILTERED/545_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

# patient 564 not hyperdiploid!
#perl $SCRIPT 564 rem_dia $VCF_DIR/564_rem_dia_calls_snpeff.vcf 564_rem 564_dia snp --vcf-out $VCF_OUT_FILTERED/564_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
#perl $SCRIPT 564 rem_dia $VCF_DIR_INDEL/564_rem_dia_snpeff.vcf 564_rem 564_dia indel --vcf-out $VCF_OUT_FILTERED/564_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
#perl $SCRIPT 564 rem_rel $VCF_DIR/564_rem_rel_calls_snpeff.vcf 564_rem 564_rel snp --vcf-out $VCF_OUT_FILTERED/564_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
#perl $SCRIPT 564 rem_rel $VCF_DIR_INDEL/564_rem_rel_snpeff.vcf 564_rem 564_rel indel --vcf-out $VCF_OUT_FILTERED/564_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 592 rem_dia $VCF_DIR/592_rem_dia_calls_snpeff.vcf 592_rem 592_dia snp --vcf-out $VCF_OUT_FILTERED/592_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 592 rem_dia $VCF_DIR_INDEL/592_rem_dia_snpeff.vcf 592_rem 592_dia indel --vcf-out $VCF_OUT_FILTERED/592_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 592 rem_rel $VCF_DIR/592_rem_rel_calls_snpeff.vcf 592_rem 592_rel snp --vcf-out $VCF_OUT_FILTERED/592_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 592 rem_rel $VCF_DIR_INDEL/592_rem_rel_snpeff.vcf 592_rem 592_rel indel --vcf-out $VCF_OUT_FILTERED/592_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 715 rem_dia $VCF_DIR/715_rem_dia_calls_snpeff.vcf 715_rem 715_dia snp --vcf-out $VCF_OUT_FILTERED/715_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 715 rem_dia $VCF_DIR_INDEL/715_rem_dia_snpeff.vcf 715_rem 715_dia indel --vcf-out $VCF_OUT_FILTERED/715_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 715 rem_rel $VCF_DIR/715_rem_rel_calls_snpeff.vcf 715_rem 715_rel snp --vcf-out $VCF_OUT_FILTERED/715_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 715 rem_rel $VCF_DIR_INDEL/715_rem_rel_snpeff.vcf 715_rem 715_rel indel --vcf-out $VCF_OUT_FILTERED/715_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 786 rem_dia $VCF_DIR/786_rem_dia_calls_snpeff.vcf 786_rem 786_dia snp --vcf-out $VCF_OUT_FILTERED/786_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 786 rem_dia $VCF_DIR_INDEL/786_rem_dia_snpeff.vcf 786_rem 786_dia indel --vcf-out $VCF_OUT_FILTERED/786_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 786 rem_rel $VCF_DIR/786_rem_rel_calls_snpeff.vcf 786_rem 786_rel snp --vcf-out $VCF_OUT_FILTERED/786_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 786 rem_rel $VCF_DIR_INDEL/786_rem_rel_snpeff.vcf 786_rem 786_rel indel --vcf-out $VCF_OUT_FILTERED/786_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 792 rem_dia $VCF_DIR/792_rem_dia_calls_snpeff.vcf 792_rem 792_dia snp --vcf-out $VCF_OUT_FILTERED/792_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 792 rem_dia $VCF_DIR_INDEL/792_rem_dia_snpeff.vcf 792_rem 792_dia indel --vcf-out $VCF_OUT_FILTERED/792_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 792 rem_rel $VCF_DIR/792_rem_rel_calls_snpeff.vcf 792_rem 792_rel snp --vcf-out $VCF_OUT_FILTERED/792_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 792 rem_rel $VCF_DIR_INDEL/792_rem_rel_snpeff.vcf 792_rem 792_rel indel --vcf-out $VCF_OUT_FILTERED/792_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 818 rem_dia $VCF_DIR/818_rem_dia_calls_snpeff.vcf 818_rem 818_dia snp --vcf-out $VCF_OUT_FILTERED/818_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 818 rem_dia $VCF_DIR_INDEL/818_rem_dia_snpeff.vcf 818_rem 818_dia indel --vcf-out $VCF_OUT_FILTERED/818_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 818 rem_rel $VCF_DIR/818_rem_rel_calls_snpeff.vcf 818_rem 818_rel snp --vcf-out $VCF_OUT_FILTERED/818_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 818 rem_rel $VCF_DIR_INDEL/818_rem_rel_snpeff.vcf 818_rem 818_rel indel --vcf-out $VCF_OUT_FILTERED/818_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 842 rem_dia $VCF_DIR/842_rem_dia_calls_snpeff.vcf 842_rem 842_dia snp --vcf-out $VCF_OUT_FILTERED/842_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 842 rem_dia $VCF_DIR_INDEL/842_rem_dia_snpeff.vcf 842_rem 842_dia indel --vcf-out $VCF_OUT_FILTERED/842_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 842 rem_rel $VCF_DIR/842_rem_rel_calls_snpeff.vcf 842_rem 842_rel snp --vcf-out $VCF_OUT_FILTERED/842_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 842 rem_rel $VCF_DIR_INDEL/842_rem_rel_snpeff.vcf 842_rem 842_rel indel --vcf-out $VCF_OUT_FILTERED/842_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 1021247 rem_dia $VCF_DIR/1021247_rem_dia_calls_snpeff.vcf 1021247_rem 1021247_dia snp --vcf-out $VCF_OUT_FILTERED/1021247_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 1021247 rem_dia $VCF_DIR_INDEL/1021247_rem_dia_snpeff.vcf 1021247_rem 1021247_dia indel --vcf-out $VCF_OUT_FILTERED/1021247_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 1021247 rem_rel $VCF_DIR/1021247_rem_rel_calls_snpeff.vcf 1021247_rem 1021247_rel snp --vcf-out $VCF_OUT_FILTERED/1021247_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 1021247 rem_rel $VCF_DIR_INDEL/1021247_rem_rel_snpeff.vcf 1021247_rem 1021247_rel indel --vcf-out $VCF_OUT_FILTERED/1021247_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT A rem_dia $VCF_DIR/A_rem_dia_calls_snpeff.vcf A13324_rem A12642_dia snp --vcf-out $VCF_OUT_FILTERED/A_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT A rem_dia $VCF_DIR_INDEL/A_rem_dia_snpeff.vcf A13324_rem A12642_dia indel --vcf-out $VCF_OUT_FILTERED/A_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT A rem_rel $VCF_DIR/A_rem_rel_calls_snpeff.vcf A13324_rem A12886_rel snp --vcf-out $VCF_OUT_FILTERED/A_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT A rem_rel $VCF_DIR_INDEL/A_rem_rel_snpeff.vcf A13324_rem A12886_rel indel --vcf-out $VCF_OUT_FILTERED/A_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT B rem_dia $VCF_DIR/B_rem_dia_calls_snpeff.vcf B20946_rem B19668_dia snp --vcf-out $VCF_OUT_FILTERED/B_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT B rem_dia $VCF_DIR_INDEL/B_rem_dia_snpeff.vcf B20946_rem B19668_dia indel --vcf-out $VCF_OUT_FILTERED/B_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT B rem_rel $VCF_DIR/B_rem_rel_calls_snpeff.vcf B20946_rem B15010_rel snp --vcf-out $VCF_OUT_FILTERED/B_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT B rem_rel $VCF_DIR_INDEL/B_rem_rel_snpeff.vcf B20946_rem B15010_rel indel --vcf-out $VCF_OUT_FILTERED/B_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT C rem_dia $VCF_DIR/C_rem_dia_calls_snpeff.vcf C20499_rem C19797_dia snp --vcf-out $VCF_OUT_FILTERED/C_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT C rem_dia $VCF_DIR_INDEL/C_rem_dia_snpeff.vcf C20499_rem C19797_dia indel --vcf-out $VCF_OUT_FILTERED/C_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT C rem_rel $VCF_DIR/C_rem_rel_calls_snpeff.vcf C20499_rem C15050_rel snp --vcf-out $VCF_OUT_FILTERED/C_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT C rem_rel $VCF_DIR_INDEL/C_rem_rel_snpeff.vcf C20499_rem C15050_rel indel --vcf-out $VCF_OUT_FILTERED/C_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT D rem_dia $VCF_DIR/D_rem_dia_calls_snpeff.vcf D4502_rem D3826_dia snp --vcf-out $VCF_OUT_FILTERED/D_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT D rem_dia $VCF_DIR_INDEL/D_rem_dia_snpeff.vcf D4502_rem D3826_dia indel --vcf-out $VCF_OUT_FILTERED/D_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT D rem_rel $VCF_DIR/D_rem_rel_calls_snpeff.vcf D4502_rem D10183_rel snp --vcf-out $VCF_OUT_FILTERED/D_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT D rem_rel $VCF_DIR_INDEL/D_rem_rel_snpeff.vcf D4502_rem D10183_rel indel --vcf-out $VCF_OUT_FILTERED/D_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT E rem_dia $VCF_DIR/E_rem_dia_calls_snpeff.vcf E13861_rem E13174_dia snp --vcf-out $VCF_OUT_FILTERED/E_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT E rem_dia $VCF_DIR_INDEL/E_rem_dia_snpeff.vcf E13861_rem E13174_dia indel --vcf-out $VCF_OUT_FILTERED/E_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT E rem_rel $VCF_DIR/E_rem_rel_calls_snpeff.vcf E13861_rem E13479_rel snp --vcf-out $VCF_OUT_FILTERED/E_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT E rem_rel $VCF_DIR_INDEL/E_rem_rel_snpeff.vcf E13861_rem E13479_rel indel --vcf-out $VCF_OUT_FILTERED/E_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT X rem_dia $VCF_DIR/X_rem_dia_calls_snpeff.vcf X1847_rem X1286_dia snp --vcf-out $VCF_OUT_FILTERED/X_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT X rem_dia $VCF_DIR_INDEL/X_rem_dia_snpeff.vcf X1847_rem X1286_dia indel --vcf-out $VCF_OUT_FILTERED/X_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT X rem_rel $VCF_DIR/X_rem_rel_calls_snpeff.vcf X1847_rem X12831_rel snp --vcf-out $VCF_OUT_FILTERED/X_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT X rem_rel $VCF_DIR_INDEL/X_rem_rel_snpeff.vcf X1847_rem X12831_rel indel --vcf-out $VCF_OUT_FILTERED/X_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT Y rem_dia $VCF_DIR/Y_rem_dia_calls_snpeff.vcf Y3767_rem Y3141_dia snp --vcf-out $VCF_OUT_FILTERED/Y_rem_dia.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT Y rem_dia $VCF_DIR_INDEL/Y_rem_dia_snpeff.vcf Y3767_rem Y3141_dia indel --vcf-out $VCF_OUT_FILTERED/Y_rem_dia.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT Y rem_rel $VCF_DIR/Y_rem_rel_calls_snpeff.vcf Y3767_rem Y10284_rel snp --vcf-out $VCF_OUT_FILTERED/Y_rem_rel.snp.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT Y rem_rel $VCF_DIR_INDEL/Y_rem_rel_snpeff.vcf Y3767_rem Y10284_rel indel --vcf-out $VCF_OUT_FILTERED/Y_rem_rel.indel.filtered.vcf >>$RESULT_FILE 2>>$LOG_FILE

exit

# DAVID pathway enrichment analysis
cat ~/hdall/results/gene-patient-matrix.tsv | perl ~/git/hdall/get-smg.pl dia 1 | perl ~/git/hdall/pathway-analysis/pathway-enrichment-david.pl > ~/hdall/results/enriched-pathways-dia-minfreq1.david.tsv 2>~/hdall/results/enriched-pathways-dia-minfreq1.david.log
cat ~/hdall/results/gene-patient-matrix.tsv | perl ~/git/hdall/get-smg.pl dia 2 | perl ~/git/hdall/pathway-analysis/pathway-enrichment-david.pl > ~/hdall/results/enriched-pathways-dia-minfreq2.david.tsv 2>~/hdall/results/enriched-pathways-dia-minfreq2.david.log

cat ~/hdall/results/gene-patient-matrix.tsv | perl ~/git/hdall/get-smg.pl rel 1 | perl ~/git/hdall/pathway-analysis/pathway-enrichment-david.pl > ~/hdall/results/enriched-pathways-rel-minfreq1.david.tsv 2>~/hdall/results/enriched-pathways-rel-minfreq1.david.log
cat ~/hdall/results/gene-patient-matrix.tsv | perl ~/git/hdall/get-smg.pl rel 2 | perl ~/git/hdall/pathway-analysis/pathway-enrichment-david.pl > ~/hdall/results/enriched-pathways-rel-minfreq2.david.tsv 2>~/hdall/results/enriched-pathways-rel-minfreq2.david.log

cat ~/hdall/results/gene-patient-matrix.tier1.tsv | perl ~/git/hdall/get-smg.pl dia 1 | perl ~/git/hdall/pathway-analysis/pathway-enrichment-david.pl > ~/hdall/results/enriched-pathways-dia-tier1.david.tsv 2>~/hdall/results/enriched-pathways-dia-tier1.david.log
cat ~/hdall/results/gene-patient-matrix.tier1.tsv | perl ~/git/hdall/get-smg.pl rel 1 | perl ~/git/hdall/pathway-analysis/pathway-enrichment-david.pl > ~/hdall/results/enriched-pathways-rel-tier1.david.tsv 2>~/hdall/results/enriched-pathways-rel-tier1.david.log

# pathway-patient matrix

perl ~/hdall/scripts/get-pathway-patient-matrix.pl --enriched-pathways-dia ~/hdall/results/enriched-pathways-dia-minfreq1.david.tsv --enriched-pathways-rel ~/hdall/results/enriched-pathways-rel-minfreq1.david.tsv --gene-patient-matrix ~/hdall/results/gene-patient-matrix.tsv > ~/hdall/results/pathway-patient-matrix.tsv