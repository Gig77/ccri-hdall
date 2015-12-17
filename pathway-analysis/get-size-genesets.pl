use warnings FATAL => qw( all );
use strict;

use lib "/mnt/projects/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Carp;

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

my (%canonical);
open(G,"/mnt/projects/generic/data/hg19/hg19.knownCanonical.txt") or die "could not open file /mnt/projects/generic/data/hg19/hg19.knownCanonical.txt";
while(<G>)
{
	chomp;
	my ($chrom, $chromStart, $chromEnd, $clusterId, $transcript, $protein) = split(/\t/);
	
	my $geneSymbol = $id2sym{$transcript};
	$canonical{$geneSymbol} = $transcript;
}
close(G);
INFO(scalar(keys(%canonical))." canonical genes read from file /mnt/projects/generic/data/hg19/hg19.knownCanonical.txt");

my %sym2size;
my $lines = 0;
open(G,"/mnt/projects/generic/data/hg19/hg19.knownGene.txt") or die "could not open file /mnt/projects/generic/data/hg19/hg19.knownGene.txt";
while(<G>)
{
	chomp;
	my ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $proteinID, $alignID) = split(/\t/);

	$lines++;
	my $geneSymbol = $id2sym{$name};
	
	next if (exists $canonical{$geneSymbol} and $canonical{$geneSymbol} ne $name and exists $sym2size{$geneSymbol}); # prefer canonical transcript (if available)
		
	my @es = split(",", $exonStarts);
	my @ee = split(",", $exonEnds);

	# compute cds length	
	if ($cdsStart and $cdsStart < $cdsEnd)
	{
		#print "$strand\t$cdsStart\t$cdsEnd\t$exonStarts\t$exonEnds\t";
		my ($st, $en, $cdslen);		
		for (my $i = 0; $i < @es and $cdsEnd > $es[$i]; $i ++)
		{
			next if ($cdsStart > $ee[$i]);
			$st = ($cdsStart > $es[$i] and $cdsStart < $ee[$i]) ? $cdsStart : $es[$i];
			$en = ($cdsEnd > $es[$i] and $cdsEnd < $ee[$i]) ? $cdsEnd : $ee[$i];
			$cdslen += $en-$st;
		}
	
		$sym2size{$geneSymbol} = $cdslen;
	}
}
close(G);
INFO("$lines genes read from file /mnt/projects/generic/data/hg19/hg19.knownGene.txt");

my $percentil_size = keys(%sym2size) / 100;
my @sorted = sort {$sym2size{$b} <=> $sym2size{$a}} (keys(%sym2size));

print "gene\tsize\tpercentile\n";
for (my $i = 0; $i < @sorted; $i ++)
{
	my $percentile = int($i/$percentil_size)+1;
	last if ($percentile > 1);
	print "$sorted[$i]\t", $sym2size{$sorted[$i]}, "\t", int($i/$percentil_size)+1, "\n";
}
