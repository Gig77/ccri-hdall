use warnings FATAL => qw( all );
use strict;

use lib "/mnt/projects/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Data::Dumper;
use Carp;
use Getopt::Long;

# parse detailed results first
my ($cds, $density, $format);
GetOptions
(
	"cds" => \$cds,
	"density=s" => \$density,
	"format=s" => \$format   
);

die "ERROR: invalid --format: use 'nextera' or 'haloplex'\n" if (!$format or ($format ne 'nextera' and $format ne 'haloplex'));
die "ERROR: --density not specified: 'Standard' or 'Dense'\n" if ($format eq 'nextera' and !$density);

# read id mapping
my %id2sym;
open(M, "/mnt/projects/hdall/results/id-mappings.tsv") or croak "ERROR: could not read id mappings\n";
while(<M>)
{
	chomp;
	my ($sym, $id) = split(/\t/);
	$id2sym{$id} = $sym;
}
close(M);
INFO(scalar(keys(%id2sym))." id mappgins read from file /mnt/projects/hdall/results/id-mappings.tsv");

my (%canonical, %sym2size);
open(G,"/mnt/projects/generic/data/hg19/hg19.knownCanonical.txt") or die "could not open file /mnt/projects/generic/data/hg19/hg19.knownCanonical.txt";
while(<G>)
{
	chomp;
	my ($chrom, $chromStart, $chromEnd, $clusterId, $transcript, $protein) = split(/\t/);
	
	my $geneSymbol = $id2sym{$transcript};
	my $size = $chromEnd-$chromStart;
	
	# if multiple canonical transcripts for this gene symbol, use larger one
	next if ($canonical{$geneSymbol} and $sym2size{$geneSymbol} > $size); 
	
	$canonical{$geneSymbol} = $transcript;
	$sym2size{$geneSymbol} = $size;
}
close(G);
INFO(scalar(keys(%canonical))." canonical genes read from file /mnt/projects/generic/data/hg19/hg19.knownCanonical.txt");

my $lines = 0;
my %sym2info;
open(G,"/mnt/projects/generic/data/hg19/hg19.knownGene.txt") or die "could not open file /mnt/projects/generic/data/hg19/hg19.knownGene.txt";
while(<G>)
{
	chomp;
	my ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $proteinID, $alignID) = split(/\t/);

	$lines++;
	my $geneSymbol = $id2sym{$name};
	
	next if (exists $canonical{$geneSymbol} and $canonical{$geneSymbol} ne $name); # prefer canonical transcript (if available)

	$sym2info{$geneSymbol}{chr} = "$chrom";
		
	my @es = split(",", $exonStarts);
	my @ee = split(",", $exonEnds);

	# transcript length
	{
		my $trlen;		
		for (my $i = 0; $i < @es; $i ++)
		{
			$trlen += $ee[$i]-$es[$i];
		}

		#$sym2info{$prev2sym{$geneSymbol}}{'trlen'} = $trlen if ($prev2sym{$geneSymbol});
#		$sym2info{$geneSymbol}{'trlen'} = $trlen;
	}
	
	my ($st, $en);		
	for (my $i = 0; $i < @es and $cdsEnd > $es[$i]; $i ++)
	{
		next if ($cds and $cdsStart > $ee[$i]);
			
		$st = ($cds and $cdsStart > $es[$i] and $cdsStart < $ee[$i]) ? $cdsStart : $es[$i];
		$en = ($cds and $cdsEnd > $es[$i] and $cdsEnd < $ee[$i]) ? $cdsEnd : $ee[$i];
			
		$sym2info{$geneSymbol}{exons}{$i} = "$st-$en";
	}	
}
close(G);
INFO("$lines genes read from file /mnt/projects/generic/data/hg19/hg19.knownGene.txt");

print "Chromosome,StartCoordinate,StopCoordinate,TargetType,Density,Labels\n" if ($format eq 'nextera');

my %genes;
my $written = 0;
#<STDIN>; # skip header
while(<STDIN>)
{
	chomp;
	my ($geneSymbol) = split("\t");

	if (!exists $sym2info{$geneSymbol})
	{
		WARN("Gene $geneSymbol not found");
		next;
	}
	

	foreach my $e (keys(%{$sym2info{$geneSymbol}{exons}}))
	{
		my ($st,$en) = split("-", $sym2info{$geneSymbol}{exons}{$e});
		
		if ($format eq 'nextera')
		{
			print $sym2info{$geneSymbol}{chr}, ",";
			print $st, ",";
			print $en, ",";
			print "Exon", ",";
			print "$density", ",";
			print "$geneSymbol:", $e+1, "\n";			
		}
		elsif ($format eq 'haloplex')
		{
			print $sym2info{$geneSymbol}{chr}, ":";
			print $st, "-";
			print $en, "\t";
			print "$geneSymbol";
			print " # exon ", $e+1, "\n";						
		}
	}
}


