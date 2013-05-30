use warnings;
use strict;

use Data::Dumper;
use Carp;

my $debug = 1;
$| = 1; # turn on autoflush


# read id mapping
my %id2sym;
open(M, "$ENV{HOME}/hdall/results/id-mappings.tsv") or croak "ERROR: could not read id mappings\n";
while(<M>)
{
	chomp;
	my ($sym, $id) = split(/\t/);
	$id2sym{$id} = $sym;
}
close(M);

# read gene description
my %sym2info;
open(G,"$ENV{HOME}/hdall/data/hg19/hg19.kgXref.txt") or die "could not open file $ENV{HOME}/hdall/data/hg19/hg19.kgXref.txt";
while(<G>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refSeq, $protAcc, $description, $rfamAcc, $tRnaName) = split(/\t/);

	my $sym = $id2sym{$kgID};
	$sym2info{$sym}{'description'} = $description;	
}
close(G);

my (%canonical, %sym2size);
open(G,"$ENV{HOME}/hdall/data/hg19/hg19.knownCanonical.txt") or die "could not open file $ENV{HOME}/hdall/data/hg19/hg19.knownCanonical.txt";
<G>; # skip header
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

open(G,"$ENV{HOME}/hdall/data/hg19/hg19.knownGene.txt") or die "could not open file $ENV{HOME}/hdall/data/hg19/hg19.knownGene.txt";
while(<G>)
{
	chomp;
	my ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $proteinID, $alignID) = split(/\t/);

	my $geneSymbol = $id2sym{$name};
	
	next if (exists $canonical{$geneSymbol} and $canonical{$geneSymbol} ne $name and exists $sym2info{$geneSymbol}{'cdslen'}); # prefer canonical transcript (if available)
		
	#$sym2info{$prev2sym{$geneSymbol}}{'exons'} = $exonCount if ($prev2sym{$geneSymbol});
	$sym2info{$geneSymbol}{'exons'} = $exonCount;

	my @es = split(",", $exonStarts);
	my @ee = split(",", $exonEnds);

	print STDERR "$name: $exonStarts\n" if ($geneSymbol eq "LDB3");
	print STDERR "$name: $exonEnds\n" if ($geneSymbol eq "LDB3");
	
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
	
		print STDERR "$name: $cdslen\n" if ($geneSymbol eq "LDB3");

		#$sym2info{$prev2sym{$geneSymbol}}{'cdslen'} = $cdslen if ($prev2sym{$geneSymbol});
		$sym2info{$geneSymbol}{'cdslen'} = $cdslen;
	}
}
close(G);

my %genes;
<>; # skip header
while(<>)
{
	chomp;
	my ($patient, $sample, $var_type, $chr, $pos, $dbSNP, $ref, $alt, $gene, $add_genes, $impact, $effect, 
		$dp_rem_tot, $dp_rem_ref, $dp_rem_var, $dp_leu_tot, $dp_leu_ref, $dp_leu_var, $freq, $snpeff) = split("\t");

	die "ERROR: $0: snpeff annotation missing from following line:\n$_\n"
		if (!$snpeff);

	next if ($effect =~ /(DOWNSTREAM|INTERGENIC|INTRON|UPSTREAM)/);
	
	$genes{$patient}{$sample}{$gene}{"$chr:$pos:$ref->$alt:$freq:$impact:$effect"} = $genes{$patient}{$sample}{$gene}{"$chr:$pos:$ref->$alt:$freq:$impact:$effect"} 
		? $genes{$patient}{$sample}{$gene}{"$chr:$pos:$ref->$alt:$freq:$impact:$effect"}.",$snpeff"
		: $snpeff;
}

print "patient\tcomparison\tgene\ttr_len\tcds_len\texons\tdesc\tnum_mut\tmut_effects\n";
foreach my $p (keys(%genes))
{
	foreach my $s (keys(%{$genes{$p}}))
	{
		my @sorted_genes = sort {values(%{$genes{$b}}) <=> values(%{$genes{$a}})} keys(%{$genes{$p}{$s}});
		foreach my $g (@sorted_genes)
		{
			my $info = $sym2info{$g}
				or 	print STDERR "$0: WARNING: Could not map gene $g\n";
						
			print $p, "\t", $s, "\t", $g, "\t";
			print "".($info->{'trlen'} ? $info->{'trlen'} : "")."\t";
			print "".($info->{'cdslen'} ? $info->{'cdslen'} : "")."\t";
			print "".($info->{'exons'} ? $info->{'exons'} : "")."\t";
			print "".($info->{'description'} ? $info->{'description'} : "")."\t";
			
			print scalar(values(%{$genes{$p}{$s}{$g}})), "\t";
			my $first = 1;
			foreach my $v (keys(%{$genes{$p}{$s}{$g}}))
			{
				print ";" if (!$first);
				print $v,"[",$genes{$p}{$s}{$g}{$v},"]";
				$first = 0;
			}
			print "\n";
		}		
	}
}
	
