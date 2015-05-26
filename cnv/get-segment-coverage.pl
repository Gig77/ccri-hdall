use warnings FATAL => qw( all );;
use strict;

use Getopt::Long;

my %diploid_chr = (
#	"92D" => "chr3"
);

my ($sample, $bin_size);
GetOptions
(
	"sample=s" => \$sample,
	"bin-size=i" => \$bin_size
);
die "ERROR: --bin-size not specified" if (!$bin_size);
die "ERROR: --sample not specified" if (!$sample);

print STDERR "Diploid chromosome for sample $sample: ".$diploid_chr{$sample}."\n" if (defined $diploid_chr{$sample});

# read chromosome sizes
my %chr_size;
open(IN, "/mnt/projects/generic/data/hg19/ucsc.hg19.chrom.sizes") or die "could not open ucsc.hg19.chrom.sizes\n";
while(<IN>)
{
	my ($chr, $size) = split("\t");
	next if ($chr !~ /^chr[\dXY]/);
	$chr_size{$chr} = $size;
}
close(IN);

my %bins;
my $diploid_total = 0;
while(<>)
{
    next if /^#/;
    next if /^$/;

	my ($chr, $pos, $cov) = split /\t/;
	$chr = "chr$chr" if ($chr !~ /^chr/);

	my $bin = int($pos/$bin_size);
	$bins{"$chr:$bin"} += $cov;
	$diploid_total += $cov if (defined $diploid_chr{$sample} and $chr eq $diploid_chr{$sample}); 
}

foreach my $chr (keys(%chr_size))
{
	for (my $i = 0; $i < $chr_size{$chr}; $i += $bin_size)
	{
		my $cov = $bins{"$chr:".int($i/$bin_size)};
		print "$chr\t$i\t".(defined $cov ? $cov / ($diploid_total > 0 ? $diploid_total / 1000 : 1) : 0)."\n";  
	}
}
