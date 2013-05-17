ROI_FILE=~/hdall/results/music/ucsc-genes.hg19.roi

# download exon BED file from UCSC table browser (+2 flanking bp to include splice sites)
# download human reference genome (hg19) from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/hg19/ucsc.hg19.fasta.gz

# filter and sort UCSC bed file
#grep -vP '^(\S+random|chrUn)' ~/hdall/results/ucsc-genes.hg19.bed > ~/hdall/results/ucsc-genes.hg19.filtered.bed 
~/tools/lh3-sort/sort -k 1,1N -k 2,2n ~/hdall/results/ucsc-genes.hg19.bed > ~/hdall/results/ucsc-genes.hg19.sorted.bed

# convert to 1-based start coordinate
#cat ~/hdall/results/ucsc-genes.hg19.sorted.bed | perl -ne 'if (/^([^\t]+)\t(\d+)\t(.*)/) { print "$1\t",$2+1,\t$3\n" } else { print "$_" };' > ~/hdall/results/ucsc-genes.hg19.sorted.1based.bed 

# merge overlapping exons using bedtools
~/tools/bedtools-2.17.0/bin/bedtools merge -nms -n -i ~/hdall/results/ucsc-genes.hg19.sorted.bed > ~/hdall/results/ucsc-genes.hg19.sorted.merged.bed

# postprocess ROI (map back to gene symbols, 1-based start coordinate, cut coordinates back to chromosome length)
cat ~/hdall/results/ucsc-genes.hg19.sorted.merged.bed | perl ~/git/hdall/music/ucsc-exon-to-gene-symbol.pl > $ROI_FILE
echo "Region of interest (ROI) file written to $ROI_FILE" 

# compress and index it for fast access with tabix
bgzip -c $ROI_FILE > $ROI_FILE.gz
tabix $ROI_FILE.gz -s 1 -b 2 -e 3

exit

# create WIG file compatible with MuSiC
#perl ~/git/hdall/bedgraph-to-wig.pl /home/STANNANET/christian.frech/hdall/data/bam/715_rem.merged.duplicate_marked.realigned.recalibrated.bedgraph /home/STANNANET/christian.frech/hdall/data/bam/715_rem.merged.duplicate_marked.realigned.recalibrated.music.wig
perl ~/git/hdall/bedgraph-to-wig.pl /home/STANNANET/christian.frech/hdall/data/bam/715_dia.merged.duplicate_marked.realigned.recalibrated.bedgraph /home/STANNANET/christian.frech/hdall/data/bam/715_dia.merged.duplicate_marked.realigned.recalibrated.music.wig

# run music bmr calc-wig-covg 
genome music bmr calc-wig-covg --wig-list /mnt/suse/data/christian/hdall/results/current/music/wig-list.tsv --output-dir /mnt/suse/data/christian/hdall/results/current/music/out/ --reference-sequence /mnt/suse/home/STANNANET/christian.frech/hdall/data/hg19/ucsc.hg19.fasta --roi-file /mnt/suse/data/christian/hdall/results/current/music/ucsc-genes.hg19.roi