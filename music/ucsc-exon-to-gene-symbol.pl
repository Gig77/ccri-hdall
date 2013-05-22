#!/usr/bin/perl

use strict;
use warnings;

# read translation table

#FORMAT:
#kgID	mRNA	spID	spDisplayID	geneSymbol	refseq	protAcc	description	rfamAcc	tRnaName
#uc001aaa.3	BC032353			DDX11L1			Homo sapiens mRNA for DEAD/H box polypeptide 11 like 1 (DDX11L1 gene).		
my %ucsc2sym;
open(IN,"/home/STANNANET/christian.frech/hdall/data/hg19/hg19.kgXref.txt") or die "could not open kgXref\n";
while(<IN>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refseq, $protAcc, $description, $rfamAcc, $tRnaName) = split("\t");
	$ucsc2sym{$kgID} = $geneSymbol;
}
close(IN);

# read chromosome sizes
my %chrsize;
open(IN, 
"/home/STANNANET/christian.frech/hdall/data/hg19/ucsc.hg19.chrom.sizes") or die "could not open ucsc.hg19.chrom.sizes\n";
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

	$start = $start + 1; # convert to 1-based start coordinate
	$start = $chrsize{$chr} if ($start > $chrsize{$chr}); # trim coordinate to chromosome size
	$end = $chrsize{$chr} if ($end > $chrsize{$chr}); # trim coordinate to chromosome size
	next if ($start >= $end);
	next if ($chr eq "M" and $end == 16571); # gives music error: ERROR: Request for chrM:16572-16573 in /mnt/suse/home/STANNANET/christian.frech/hdall/data/hg19/ucsc.hg19.fasta, but chrM has length 16571

	my @regions;
	while($region =~ /(uc.*?)_exon/g) { push(@regions, $1) }
	@regions = sort(@regions);
	
	my %sym_written;
	foreach my $r (@regions)
	{
		my $sym = $ucsc2sym{$r} ? $ucsc2sym{$r} : $r;
		next if ($sym_written{$sym});
		print "$chr\t$start\t$end\t$sym\t$num_regions\n"; 
		$sym_written{$sym} = 1;
	}
}
