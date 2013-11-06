use warnings FATAL => qw( all );
use strict;

# read id mapping
my %id2sym;
open(M, "$ENV{HOME}/hdall/results/id-mappings.tsv") or die "ERROR: could not read id mappings\n";
while(<M>)
{
	chomp;
	my ($sym, $id) = split(/\t/);
	$id2sym{$id} = $sym;
}
close(M);

my (%canonical, %sym2size);
open(G,"$ENV{HOME}/generic/data/hg19/hg19.knownCanonical.txt") or die "could not open file $ENV{HOME}/generic/data/hg19/hg19.knownCanonical.txt";
while(<G>)
{
	chomp;
	my ($chrom, $chromStart, $chromEnd, $clusterId, $transcript, $protein) = split(/\t/);

	my $geneSymbol = $id2sym{$transcript};
	my $size = $chromEnd-$chromStart;
	next if ($canonical{$geneSymbol} and $sym2size{$geneSymbol} > $size); # if multiple canonical transcripts for this gene symbol, use larger one 

	$canonical{$geneSymbol} = $transcript;
	$sym2size{$geneSymbol} = $size;
}
close(G);

my %sym2info;
open(G,"$ENV{HOME}/generic/data/hg19/hg19.knownGene.txt") or die "could not open file $ENV{HOME}/generic/data/hg19/hg19.knownGene.txt";
while(<G>)
{
	chomp;
	my ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $proteinID, $alignID) = split(/\t/);

	my $geneSymbol = $id2sym{$name};
	#next if (!exists $canonical{$geneSymbol} or $canonical{$geneSymbol} ne $name); # skip non-canoncial gene models
		
	$sym2info{$geneSymbol}{'exons'} = $exonCount;
	$sym2info{$geneSymbol}{'chr'} = $chrom;
	$sym2info{$geneSymbol}{'start'} = $txStart;
	$sym2info{$geneSymbol}{'end'} = $txEnd;

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
		$sym2info{$geneSymbol}{'trlen'} = $trlen;
	}
	
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
	
		#$sym2info{$prev2sym{$geneSymbol}}{'cdslen'} = $cdslen if ($prev2sym{$geneSymbol});
		$sym2info{$geneSymbol}{'cdslen'} = $cdslen;
	}
}
close(G);

print "HGNC\tcds_len\n";
foreach my $g (keys(%sym2info))
{
	print "$g", "\t", defined $sym2info{$g}{'cdslen'} ? $sym2info{$g}{'cdslen'} : "NA", "\n";
}