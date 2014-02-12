use warnings FATAL => qw( all );;
use strict;

use Getopt::Long;

my ($bin_size);
GetOptions
(
	"bin-size=i" => \$bin_size
);
die "ERROR: --bin-size not specified" if (!$bin_size);

# read chromosome sizes
my %chr_size;
open(IN, "/home/STANNANET/christian.frech/generic/data/hg19/ucsc.hg19.chrom.sizes") or die "could not open ucsc.hg19.chrom.sizes\n";
while(<IN>)
{
	my ($chr, $size) = split("\t");
	next if ($chr !~ /^chr[\dXY]/);
	$chr_size{$chr} = $size;
}
close(IN);

my %bins;
while(<>)
{
    next if /^#/;
    next if /^$/;

	my ($chr, $pos, $cov) = split /\t/;

	my $bin = int($pos/$bin_size);
	$bins{"$chr:$bin"} += $cov; 
}

foreach my $chr (keys(%chr_size))
{
	for (my $i = 0; $i < $chr_size{$chr}; $i += $bin_size)
	{
		my $cov = $bins{"$chr:".int($i/$bin_size)};
		print "$chr\t$i\t".(defined $cov ? $cov : 0)."\n";  
	}
}
