use warnings FATAL => qw( all );
use strict;
use Carp;

use lib "/mnt/projects/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Vcf;
use Getopt::Long;
use Data::Dumper;

#Log::Log4perl->get_logger()->level($ERROR);

INFO("START: ", join(" ", @ARGV));

#my ($music_roi, $keep_variants_file, $entrez_mapping, $sample, $header, $min_af, $deleterious);
#GetOptions
#(
#	"music-roi=s" => \$music_roi,  # MuSiC region of interest (ROI) file (must be tabix accessible, i.e. compressed and indexed)
##	"keep-variants-file=s" => \$keep_variants_file,  # tab-separated file with variants to keep (chr, start)
#	"mapping-entrez=s" => \$entrez_mapping,  # file with mappings from gene symbol to entrez ids
#	"sample=s" => \$sample,  # e.g. 314_rem_dia
#	"header" => \$header,  # output header yes/no
#	"deleterious" => \$deleterious,  # output only variants predicted to be deleterious by PolyPhen or SIFT
#	"min-af=s" => \$min_af  # minimum allelic frequency
#);


my $vcf = Vcf->new(file => "-");
$vcf->parse_header();
my (@samples) = $vcf->get_samples();

while (my $x = $vcf->next_data_hash())
{	
	my $chr = $x->{CHROM};
	my $pos = $x->{POS};
	my $ref_allele = $x->{REF};

	##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
	##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as d

	foreach my $sample (@samples)
	{
		next if ($sample !~ /_rem/);
		
		my $gt = $x->{gtypes}{$sample}{GT};
		my $dp = $x->{gtypes}{$sample}{DP};
		my $gq = $x->{gtypes}{$sample}{GQ};
	
		next if ($gt eq "./.");	
		next if ($dp < 3);
	
		for(my $i = 0; $i < @{$x->{ALT}}; $i ++)
		{
			my $ad = (split(",", $x->{gtypes}{$sample}{AD}))[$i+1];
			next if ($ad < 2);
			
			my $alt_allele = $x->{ALT}->[$i];
			
			print "$sample\t$chr\t$pos\t$ref_allele\t$alt_allele\t$dp\t$ad\t$gt\n";
		}			
	}
}
$vcf->close();

INFO("END");
