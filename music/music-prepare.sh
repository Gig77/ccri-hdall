exit

# get (scrape) DAVID pathway annotations
# ATTENTION: I ran this script two times with different sorting of input genes, because for some reason the 
# david web service FAILS TO ANNOTATE THE FIRST GENE in each query!
# the final list is produced by merging the results from both runs 
perl /mnt/projects/hdall/scripts/get-david-annotations.pl > /mnt/projects/hdall/results/music/pathways-david.tsv
perl /mnt/projects/hdall/scripts/get-david-annotations.pl > /mnt/projects/hdall/results/music/pathways-david-reverse.tsv # change sort in script before running it!
perl /mnt/projects/hdall/scripts/merge-david-annotations.pl /mnt/projects/hdall/results/music/pathways-david.tsv /mnt/projects/hdall/results/music/pathways-david-reverse.tsv > /mnt/projects/hdall/results/music/pathways-david-complete.tsv

# get example ROI file (Ensembl v67)
cd /mnt/projects/hdall/results/music/
curl -O https://raw.github.com/ding-lab/calc-roi-covg/master/data/ensembl_67_cds_ncrna_and_splice_sites_hg19
mv ensembl_67_cds_ncrna_and_splice_sites_hg19 ensembl_67_cds_ncrna_and_splice_sites_hg19.roi
~/tools/tabix-0.2.6/bgzip /mnt/projects/hdall/results/music/ensembl_67_cds_ncrna_and_splice_sites_hg19.roi
~/tools/tabix-0.2.6/tabix /mnt/projects/hdall/results/music/ensembl_67_cds_ncrna_and_splice_sites_hg19.roi.gz -b 2 -e 3

# download exon BED file from UCSC table browser (+2 flanking bp to include splice sites)
# download human reference genome (hg19) from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/hg19/ucsc.hg19.fasta.gz

# filter and sort UCSC bed file
~/tools/lh3-sort/sort -k 1,1N -k 2,2n /mnt/projects/hdall/results/ucsc-genes.hg19.bed > /mnt/projects/hdall/results/ucsc-genes.hg19.sorted.bed

# merge overlapping exons using bedtools
/data_synology/software/bedtools-2.17.0/bin/bedtools merge -s -nms -n -i /mnt/projects/hdall/results/ucsc-genes.hg19.sorted.bed | ~/tools/lh3-sort/sort -k 1,1N -k 2,2n > /mnt/projects/hdall/results/ucsc-genes.hg19.sorted.merged.bed

# postprocess ROI (map back to gene symbols, 1-based start coordinate, cut back coordinates exceeding chromosome length)
cat /mnt/projects/hdall/results/ucsc-genes.hg19.sorted.merged.bed | perl ~/git/hdall/music/ucsc-exon-to-gene-symbol.pl > /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi
echo "Region of interest (ROI) file written to /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi" 

# compress and index it for fast access with tabix
bgzip -c /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi > /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz
tabix /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz -s 1 -b 2 -e 3

exit

# create WIG file compatible with MuSiC
#perl ~/git/hdall/get-music-compatible-wig.pl

# convert VCF to MAF

cat /mnt/projects/hdall/results/filtered_variants/314_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 314_dia --sample-normal 314_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv --header > /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/399_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 399_dia --sample-normal 399_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/430_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 430_dia --sample-normal 430_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/446_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 446_dia --sample-normal 446_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/460_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 460_dia --sample-normal 460_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/545_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 545_dia --sample-normal 545_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/592_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 592_dia --sample-normal 592_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/715_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 715_dia --sample-normal 715_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/786_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 786_dia --sample-normal 786_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/792_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 792_dia --sample-normal 792_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/818_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 818_dia --sample-normal 818_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/842_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 842_dia --sample-normal 842_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/1021247_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor 1021247_dia --sample-normal 1021247_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/A_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor A_dia --sample-normal A_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/B_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor B_dia --sample-normal B_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/C_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor C_dia --sample-normal C_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/D_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor D_dia --sample-normal D_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/E_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor E_dia --sample-normal E_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/X_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor X_dia --sample-normal X_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log
cat /mnt/projects/hdall/results/filtered_variants/Y_rem_dia.snp.filtered.vcf | perl /mnt/projects/hdall/scripts/vcf2maf.pl --sample-tumor Y_dia --sample-normal Y_rem --music-roi /mnt/projects/hdall/results/music/ucsc-genes.hg19.roi.gz --mapping-entrez /mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv >> /mnt/projects/hdall/results/current/music/dia/rem_dia.maf 2>>/mnt/projects/hdall/results/music/dia/rem_dia.maf.log

# run music bmr calc-wig-covg
# needs to be run locally, because music is currently not installed on suse 
genome music bmr calc-wig-covg --wig-list /mnt/suse/data/christian/hdall/results/current/music/dia/wig-list.tsv --output-dir /mnt/suse/data/christian/hdall/results/current/music/dia/ --reference-sequence /mnt/suse/data/christian/hdall/data/current/hg19/ucsc.hg19.nochr.fasta --roi-file /mnt/suse/data/christian/hdall/results/current/music/ucsc-genes.hg19.roi
genome music bmr calc-wig-covg --wig-list /mnt/suse/data/christian/hdall/results/current/music/rel/wig-list.tsv --output-dir /mnt/suse/data/christian/hdall/results/current/music/rel/ --reference-sequence /mnt/suse/data/christian/hdall/data/current/hg19/ucsc.hg19.nochr.fasta --roi-file /mnt/suse/data/christian/hdall/results/current/music/ucsc-genes.hg19.roi

# run music bmr calc-bmr
genome music bmr calc-bmr --bam-list /mnt/suse/data/christian/hdall/results/current/music/dia/wig-list.tsv --maf-file /mnt/suse/data/christian/hdall/results/current/music/dia/rem_dia.maf --output-dir /mnt/suse/data/christian/hdall/results/current/music/dia/ --reference-sequence /mnt/suse/data/christian/hdall/data/current/hg19/ucsc.hg19.nochr.fasta --roi-file /mnt/suse/data/christian/hdall/results/current/music/ucsc-genes.hg19.roi

# run music smg
genome music smg --gene-mr-file /mnt/suse/data/christian/hdall/results/current/music/dia/gene_mrs --output-file /mnt/suse/data/christian/hdall/results/current/music/dia/smg.tsv

# run music path-scan
genome music path-scan --bam-list /mnt/suse/data/christian/hdall/results/current/music/dia/wig-list.tsv --gene-covg-dir /mnt/suse/data/christian/hdall/results/current/music/dia/gene_covgs/ --maf-file /mnt/suse/data/christian/hdall/results/current/music/dia/rem_dia.maf --output-file /mnt/suse/data/christian/hdall/results/current/music/dia/sm_pathways.tsv --pathway-file /mnt/suse/data/christian/hdall/results/current/music/pathways-david-complete.tsv
perl /mnt/projects/hdall/scripts/pathway-analysis/annotate-pathscan-result.pl --sm-pathways /mnt/projects/hdall/results/music/dia/sm_pathways.tsv --sm-pathways-detail /mnt/projects/hdall/results/music/dia/sm_pathways.tsv_detailed > /mnt/projects/hdall/results/music/dia/sm_pathways.annotated.tsv