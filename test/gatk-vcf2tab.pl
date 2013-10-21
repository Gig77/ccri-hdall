use warnings FATAL => qw( all );
use strict;
use Carp;

use lib "$ENV{HOME}/generic/scripts";
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
	my $alt_allele = join(",", @{$x->{ALT}});
	my $dbSNP = ($x->{ID} and $x->{ID} ne '.') ? $x->{ID} : "";
	
	##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
	##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as d
	my $ad = $x->{gtypes}{'1021247_rem'}{AD};
	my $dp = $x->{gtypes}{'1021247_rem'}{DP};
	my $gq = $x->{gtypes}{'1021247_rem'}{GQ};
	my $gt = $x->{gtypes}{'1021247_rem'}{GT};
	my $pl = $x->{gtypes}{'1021247_rem'}{PL};

	next if ($gt eq "./.");
	
	my $ad_ref = (split(",", $ad))[0];
	my $ad_alt = $gt =~ /2/ ? (split(",", $ad))[2] : (split(",", $ad))[1];
	 
	next if ($dp < 3);
	
	if (!$ad)
	{
		print "$chr:$pos\n";
		print Dumper($x->{gtypes}{'1021247_rem'});
		exit;
	}
	print "$chr\t$pos\t$ref_allele\t$alt_allele\t$ad\t$dp\t$gq\t$gt\t$pl\t$ad_ref\t$ad_alt\n";
}
$vcf->close();

INFO("END");
