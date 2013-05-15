#!/usr/bin/perl

use strict;
use warnings;

# read translation table

#FORMAT:
#kgID	mRNA	spID	spDisplayID	geneSymbol	refseq	protAcc	description	rfamAcc	tRnaName
#uc001aaa.3	BC032353			DDX11L1			Homo sapiens mRNA for DEAD/H box polypeptide 11 like 1 (DDX11L1 gene).		
my %ucsc2sym;
open(IN,"/data/christian/kamilla/data/current/hg19/hg19.kgXref.txt") or die "could not open kgXref\n";
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
"/data/christian/kamilla/data/current/hg19/ucsc.hg19.chrom.sizes") or die "could not open ucsc.hg19.chrom.sizes\n";
while(<IN>)
{
	my ($chr, $size) = split("\t");
	$chrsize{$chr} = $size;
}
close(IN);

#chr1	136950	139698	uc001aam.4_exon_2_2_chr1_136953_r;uc021oeg.1_exon_0_2_chr1_137841_r	2

# first iteration: determine number of regions for each symbol
my (%symbols, @lines);
while(<>)
{
	chomp;
	push(@lines, $_); # remember for second iteration
	my ($chr, $start, $end, $region) = split("\t");
	while($region =~ /(uc.*?)_exon/g)
	{
		my $sym = $ucsc2sym{$1} ? $ucsc2sym{$1} : $1;
		$symbols{$sym} = exists $symbols{$sym} ? $symbols{$sym} + 1 : 1;
	}
}

# second iteration: select symbol with most regions as common identifier in output ROI
foreach my $line (@lines)
{
	my ($chr, $start, $end, $region) = split("\t", $line);

	$start = $start + 1; # convert to 1-based start coordinate
	$start = $chrsize{$chr} if ($start > $chrsize{$chr}); # trim coordinate to chromosome size
	$end = $chrsize{$chr} if ($end > $chrsize{$chr}); # trim coordinate to chromosome size
	next if ($start >= $end);
	next if ($chr eq "chrM" and $end == 16571); # gives music error: ERROR: Request for chrM:16572-16573 in /mnt/suse/data/christian/kamilla/data/current/hg19/ucsc.hg19.fasta, but chrM has length 16571
	
	print $chr, "\t", $start, "\t", $end, "\t"; 

	# determine which of the symbols has most target regions and use this one	
	my ($pref_symbol, $max_regions) = (undef, 0);
	my %region_symbols;
	while($region =~ /(uc.*?)_exon/g)
	{
		my $sym = $ucsc2sym{$1} ? $ucsc2sym{$1} : $1;
		$region_symbols{$sym} = 1;
		if ($symbols{$sym} > $max_regions)
		{
			$max_regions = $symbols{$sym};
			$pref_symbol = $sym;
		}
	}

	print "$pref_symbol\t",join(";", keys(%region_symbols)),"\n";
}