#!/usr/bin/perl

use strict;
use warnings FATAL => qw( all );

use Carp;
 
# read translation table (id mappings)
my %sym;
open(M, "/mnt/projects/hdall/results/id-mappings.tsv") or croak "ERROR: could not read id mappings\n";
while(<M>)
{
	chomp;
	my ($symbol, $id) = split(/\t/);
	$sym{$id} = $symbol;
}
close(M);

# read chromosome sizes
my %chrsize;
open(IN, 
"/mnt/projects/generic/data/hg19/ucsc.hg19.chrom.sizes") or die "could not open ucsc.hg19.chrom.sizes\n";
while(<IN>)
{
	my ($chr, $size) = split("\t");
    $chr =~ s/^chr//;
	$chrsize{$chr} = $size;
}
close(IN);

#chr1	136950	139698	uc001aam.4_exon_2_2_chr1_136953_r;uc021oeg.1_exon_0_2_chr1_137841_r	2
while(<>)
{
	chomp;

	my ($chr, $start, $end, $region, $num_regions) = split("\t");
    $chr =~ s/^chr//;

	$start = $start - 1 if ($start > 0); # convert to 1-based start coordinate, add padding left for splice site mutations
	$end = $end + 2; # add padding right for splice site mutations
	$start = $chrsize{$chr} if ($start > $chrsize{$chr}); # trim coordinate to chromosome size
	$end = $chrsize{$chr} if ($end > $chrsize{$chr}); # trim coordinate to chromosome size
	next if ($start >= $end);
	next if ($chr eq "M" and $end == 16571); # gives music error: ERROR: Request for chrM:16572-16573 in /mnt/suse/mnt/projects/generic/data/hg19/ucsc.hg19.fasta, but chrM has length 16571

	my @regions;
	while($region =~ /(uc.*?)_exon/g) { push(@regions, $1) }
	@regions = sort(@regions);
	
	my %sym_written;
	foreach my $r (@regions)
	{
		my $symbol = $sym{$r} ? $sym{$r} : $r;
		next if ($sym_written{$symbol});
		print "$chr\t$start\t$end\t$symbol\t$num_regions\n"; 
		$sym_written{$symbol} = 1;
	}
}
