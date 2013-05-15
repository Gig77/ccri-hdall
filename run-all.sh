SCRIPT=~/kamilla/scripts/filter-variants.pl
VCF_DIR=/home/STANNANET/christian.frech/kamilla/data/current/mutect_vcf
VCF_DIR_INDEL=/home/STANNANET/christian.frech/kamilla/data/current/somatic_indel_vcf
RESULT_FILE=~/kamilla/results/current/filtered-variants.tsv
LOG_FILE=~/kamilla/results/current/filtered-variants.log

perl $SCRIPT 314 rem_dia $VCF_DIR/314_rem_dia_calls_snpeff.vcf 314_rem 314_dia snp 1 >$RESULT_FILE 2>$LOG_FILE
#perl $SCRIPT 314 rem_dia $VCF_DIR_INDEL/314_rem_dia_snpeff.vcf 314_rem 314_dia indel >>$RESULT_FILE 2>$LOG_FILE
perl $SCRIPT 314 rem_rel $VCF_DIR/314_rem_rel_calls_snpeff.vcf 314_rem 314_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 314 rem_rel $VCF_DIR_INDEL/314_rem_rel_snpeff.vcf 314_rem 314_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 399 rem_dia $VCF_DIR/399_rem_dia_calls_snpeff.vcf 399_rem 399_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 399 rem_dia $VCF_DIR_INDEL/399_rem_dia_snpeff.vcf 399_rem 399_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 399 rem_rel $VCF_DIR/399_rem_rel_calls_snpeff.vcf 399_rem 399_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 399 rem_rel $VCF_DIR_INDEL/399_rem_rel_snpeff.vcf 399_rem 399_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 430 rem_dia $VCF_DIR/430_rem_dia_calls_snpeff.vcf 430_rem 430_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 430 rem_dia $VCF_DIR_INDEL/430_rem_dia_snpeff.vcf 430_rem 430_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 430 rem_rel $VCF_DIR/430_rem_rel_calls_snpeff.vcf 430_rem 430_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 430 rem_rel $VCF_DIR_INDEL/430_rem_rel_snpeff.vcf 430_rem 430_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 446 rem_dia $VCF_DIR/446_rem_dia_calls_snpeff.vcf 446_rem 446_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 446 rem_dia $VCF_DIR_INDEL/446_rem_dia_snpeff.vcf 446_rem 446_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 446 rem_rel $VCF_DIR/446_rem_rel_calls_snpeff.vcf 446_rem 446_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 446 rem_rel $VCF_DIR_INDEL/446_rem_rel_snpeff.vcf 446_rem 446_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 460 rem_dia $VCF_DIR/460_rem_dia_calls_snpeff.vcf 460_rem 460_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 460 rem_dia $VCF_DIR_INDEL/460_rem_dia_snpeff.vcf 460_rem 460_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 460 rem_rel $VCF_DIR/460_rem_rel_calls_snpeff.vcf 460_rem 460_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 460 rem_rel $VCF_DIR_INDEL/460_rem_rel_snpeff.vcf 460_rem 460_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 545 rem_dia $VCF_DIR/545_rem_dia_calls_snpeff.vcf 545_rem 545_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 545 rem_dia $VCF_DIR_INDEL/545_rem_dia_snpeff.vcf 545_rem 545_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 545 rem_rel $VCF_DIR/545_rem_rel_calls_snpeff.vcf 545_rem 545_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 545 rem_rel $VCF_DIR_INDEL/545_rem_rel_snpeff.vcf 545_rem 545_rel indel >>$RESULT_FILE 2>>$LOG_FILE

# patient 564 not hyperdiploid!
#perl $SCRIPT 564 rem_dia $VCF_DIR/564_rem_dia_calls_snpeff.vcf 564_rem 564_dia snp >>$RESULT_FILE 2>>$LOG_FILE
#perl $SCRIPT 564 rem_dia $VCF_DIR_INDEL/564_rem_dia_snpeff.vcf 564_rem 564_dia indel >>$RESULT_FILE 2>>$LOG_FILE
#perl $SCRIPT 564 rem_rel $VCF_DIR/564_rem_rel_calls_snpeff.vcf 564_rem 564_rel snp >>$RESULT_FILE 2>>$LOG_FILE
#perl $SCRIPT 564 rem_rel $VCF_DIR_INDEL/564_rem_rel_snpeff.vcf 564_rem 564_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 592 rem_dia $VCF_DIR/592_rem_dia_calls_snpeff.vcf 592_rem 592_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 592 rem_dia $VCF_DIR_INDEL/592_rem_dia_snpeff.vcf 592_rem 592_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 592 rem_rel $VCF_DIR/592_rem_rel_calls_snpeff.vcf 592_rem 592_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 592 rem_rel $VCF_DIR_INDEL/592_rem_rel_snpeff.vcf 592_rem 592_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 715 rem_dia $VCF_DIR/715_rem_dia_calls_snpeff.vcf 715_rem 715_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 715 rem_dia $VCF_DIR_INDEL/715_rem_dia_snpeff.vcf 715_rem 715_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 715 rem_rel $VCF_DIR/715_rem_rel_calls_snpeff.vcf 715_rem 715_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 715 rem_rel $VCF_DIR_INDEL/715_rem_rel_snpeff.vcf 715_rem 715_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 786 rem_dia $VCF_DIR/786_rem_dia_calls_snpeff.vcf 786_rem 786_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 786 rem_dia $VCF_DIR_INDEL/786_rem_dia_snpeff.vcf 786_rem 786_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 786 rem_rel $VCF_DIR/786_rem_rel_calls_snpeff.vcf 786_rem 786_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 786 rem_rel $VCF_DIR_INDEL/786_rem_rel_snpeff.vcf 786_rem 786_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 792 rem_dia $VCF_DIR/792_rem_dia_calls_snpeff.vcf 792_rem 792_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 792 rem_dia $VCF_DIR_INDEL/792_rem_dia_snpeff.vcf 792_rem 792_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 792 rem_rel $VCF_DIR/792_rem_rel_calls_snpeff.vcf 792_rem 792_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 792 rem_rel $VCF_DIR_INDEL/792_rem_rel_snpeff.vcf 792_rem 792_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 818 rem_dia $VCF_DIR/818_rem_dia_calls_snpeff.vcf 818_rem 818_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 818 rem_dia $VCF_DIR_INDEL/818_rem_dia_snpeff.vcf 818_rem 818_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 818 rem_rel $VCF_DIR/818_rem_rel_calls_snpeff.vcf 818_rem 818_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 818 rem_rel $VCF_DIR_INDEL/818_rem_rel_snpeff.vcf 818_rem 818_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 842 rem_dia $VCF_DIR/842_rem_dia_calls_snpeff.vcf 842_rem 842_dia snp >>$RESULT_FILE 2>>$LOG_FILE
#perl $SCRIPT 842 rem_dia $VCF_DIR_INDEL/842_rem_dia_snpeff.vcf 842_rem 842_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 842 rem_rel $VCF_DIR/842_rem_rel_calls_snpeff.vcf 842_rem 842_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 842 rem_rel $VCF_DIR_INDEL/842_rem_rel_snpeff.vcf 842_rem 842_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT 1021247 rem_dia $VCF_DIR/1021247_rem_dia_calls_snpeff.vcf 1021247_rem 1021247_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 1021247 rem_dia $VCF_DIR_INDEL/1021247_rem_dia_snpeff.vcf 1021247_rem 1021247_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 1021247 rem_rel $VCF_DIR/1021247_rem_rel_calls_snpeff.vcf 1021247_rem 1021247_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT 1021247 rem_rel $VCF_DIR_INDEL/1021247_rem_rel_snpeff.vcf 1021247_rem 1021247_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT A rem_dia $VCF_DIR/A_rem_dia_calls_snpeff.vcf A13324_rem A12642_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT A rem_dia $VCF_DIR_INDEL/A_rem_dia_snpeff.vcf A13324_rem A12642_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT A rem_rel $VCF_DIR/A_rem_rel_calls_snpeff.vcf A13324_rem A12886_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT A rem_rel $VCF_DIR_INDEL/A_rem_rel_snpeff.vcf A13324_rem A12886_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT B rem_dia $VCF_DIR/B_rem_dia_calls_snpeff.vcf B20946_rem B19668_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT B rem_dia $VCF_DIR_INDEL/B_rem_dia_snpeff.vcf B20946_rem B19668_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT B rem_rel $VCF_DIR/B_rem_rel_calls_snpeff.vcf B20946_rem B15010_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT B rem_rel $VCF_DIR_INDEL/B_rem_rel_snpeff.vcf B20946_rem B15010_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT C rem_dia $VCF_DIR/C_rem_dia_calls_snpeff.vcf C20499_rem C19797_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT C rem_dia $VCF_DIR_INDEL/C_rem_dia_snpeff.vcf C20499_rem C19797_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT C rem_rel $VCF_DIR/C_rem_rel_calls_snpeff.vcf C20499_rem C15050_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT C rem_rel $VCF_DIR_INDEL/C_rem_rel_snpeff.vcf C20499_rem C15050_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT D rem_dia $VCF_DIR/D_rem_dia_calls_snpeff.vcf D4502_rem D3826_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT D rem_dia $VCF_DIR_INDEL/D_rem_dia_snpeff.vcf D4502_rem D3826_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT D rem_rel $VCF_DIR/D_rem_rel_calls_snpeff.vcf D4502_rem D10183_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT D rem_rel $VCF_DIR_INDEL/D_rem_rel_snpeff.vcf D4502_rem D10183_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT E rem_dia $VCF_DIR/E_rem_dia_calls_snpeff.vcf E13861_rem E13174_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT E rem_dia $VCF_DIR_INDEL/E_rem_dia_snpeff.vcf E13861_rem E13174_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT E rem_rel $VCF_DIR/E_rem_rel_calls_snpeff.vcf E13861_rem E13479_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT E rem_rel $VCF_DIR_INDEL/E_rem_rel_snpeff.vcf E13861_rem E13479_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT X rem_dia $VCF_DIR/X_rem_dia_calls_snpeff.vcf X1847_rem X1286_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT X rem_dia $VCF_DIR_INDEL/X_rem_dia_snpeff.vcf X1847_rem X1286_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT X rem_rel $VCF_DIR/X_rem_rel_calls_snpeff.vcf X1847_rem X12831_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT X rem_rel $VCF_DIR_INDEL/X_rem_rel_snpeff.vcf X1847_rem X12831_rel indel >>$RESULT_FILE 2>>$LOG_FILE

perl $SCRIPT Y rem_dia $VCF_DIR/Y_rem_dia_calls_snpeff.vcf Y3767_rem Y3141_dia snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT Y rem_dia $VCF_DIR_INDEL/Y_rem_dia_snpeff.vcf Y3767_rem Y3141_dia indel >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT Y rem_rel $VCF_DIR/Y_rem_rel_calls_snpeff.vcf Y3767_rem Y10284_rel snp >>$RESULT_FILE 2>>$LOG_FILE
perl $SCRIPT Y rem_rel $VCF_DIR_INDEL/Y_rem_rel_snpeff.vcf Y3767_rem Y10284_rel indel >>$RESULT_FILE 2>>$LOG_FILE

# summarize variants
cat ~/kamilla/results/current/filtered-variants.tsv | perl ~/kamilla/scripts/impacted-genes.pl > ~/kamilla/results/current/impacted-genes-list.tsv
cat ~/kamilla/results/current/impacted-genes-list.tsv | perl ~/kamilla/scripts/get-gene-patient-matrix.pl > ~/kamilla/results/current/gene-patient-matrix.tsv

# pathway enrichment analysis
cat ~/kamilla/results/current/gene-patient-matrix.tsv | perl ~/kamilla/scripts/get-smg.pl dia 1 | perl ~/kamilla/scripts/pathway-enrichment-david.pl | tee ~/kamilla/results/current/enriched-pathways-dia-minfreq1.david.tsv
cat ~/kamilla/results/current/gene-patient-matrix.tsv | perl ~/kamilla/scripts/get-smg.pl dia 2 | perl ~/kamilla/scripts/pathway-enrichment-david.pl | tee ~/kamilla/results/current/enriched-pathways-dia-minfreq2.david.tsv

