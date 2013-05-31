use warnings;
use strict;

use lib "$ENV{HOME}/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
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
INFO(scalar(keys(%id2sym))." id mappgins read from file $ENV{HOME}/hdall/results/id-mappings.tsv");

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
INFO(scalar(keys(%sym2info))." gene descriptions read from file $ENV{HOME}/hdall/data/hg19/hg19.kgXref.txt");

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
INFO(scalar(keys(%canonical))." canonical genes read from file $ENV{HOME}/hdall/data/hg19/hg19.knownCanonical.txt");

my $lines = 0;
open(G,"$ENV{HOME}/hdall/data/hg19/hg19.knownGene.txt") or die "could not open file $ENV{HOME}/hdall/data/hg19/hg19.knownGene.txt";
while(<G>)
{
	chomp;
	my ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $proteinID, $alignID) = split(/\t/);

	$lines++;
	my $geneSymbol = $id2sym{$name};
	
	next if (exists $canonical{$geneSymbol} and $canonical{$geneSymbol} ne $name and exists $sym2info{$geneSymbol}{'cdslen'}); # prefer canonical transcript (if available)
		
	#$sym2info{$prev2sym{$geneSymbol}}{'exons'} = $exonCount if ($prev2sym{$geneSymbol});
	$sym2info{$geneSymbol}{'exons'} = $exonCount;

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
INFO("$lines genes read from file $ENV{HOME}/hdall/data/hg19/hg19.knownGene.txt");

my %genes;
my $written = 0;
<>; # skip header
while(<>)
{
	chomp;
	my ($patient, $sample, $var_type, $chr, $pos, $dbSNP, $ref, $alt, $gene, $add_genes, $impact, $effect, 
		$dp_rem_tot, $dp_rem_ref, $dp_rem_var, $dp_leu_tot, $dp_leu_ref, $dp_leu_var, $freq, $snpeff) = split("\t");

	die "ERROR: $0: snpeff annotation missing from following line:\n$_\n"
		if (!$snpeff);

	next if ($effect =~ /(DOWNSTREAM|INTERGENIC|INTRON|UPSTREAM|INTERGENIC_CONSERVED)/);
	
	if ($effect !~ /^(SYNONYMOUS_START|SYNONYMOUS_CODING|SYNONYMOUS_STOP|UTR_5_PRIME|UTR_5_DELETED|START_GAINED|UTR_3_PRIME|UTR_3_DELETED|INTRON_CONSERVED|INTRAGENIC|EXON)$/)
	{
		$genes{$patient}{$sample}{$gene}{'nonsyn'}{"$chr:$pos:$ref->$alt:$freq:$impact:$effect"} = $genes{$patient}{$sample}{$gene}{'nonsyn'}{"$chr:$pos:$ref->$alt:$freq:$impact:$effect"} 
			? $genes{$patient}{$sample}{$gene}{'nonsyn'}{"$chr:$pos:$ref->$alt:$freq:$impact:$effect"}.",$snpeff"
			: $snpeff;		
	}
	
	$genes{$patient}{$sample}{$gene}{'all'}{"$chr:$pos:$ref->$alt:$freq:$impact:$effect"} = $genes{$patient}{$sample}{$gene}{'all'}{"$chr:$pos:$ref->$alt:$freq:$impact:$effect"} 
		? $genes{$patient}{$sample}{$gene}{'all'}{"$chr:$pos:$ref->$alt:$freq:$impact:$effect"}.",$snpeff"
		: $snpeff;
}

print "patient\tcomparison\tgene\ttr_len\tcds_len\texons\tdesc\tnum_mut\tnum_mut_nonsyn\tmut_effects\n";
foreach my $p (keys(%genes))
{
	foreach my $s (keys(%{$genes{$p}}))
	{
		my @sorted_genes = sort {values(%{$genes{$p}{$s}{$b}{'all'}}) <=> values(%{$genes{$p}{$s}{$a}{'all'}})} keys(%{$genes{$p}{$s}});
		foreach my $g (@sorted_genes)
		{
			my $info = $sym2info{$g}
				or WARN("Could not map gene $g");
						
			print $p, "\t", $s, "\t", $g, "\t";
			print "".($info->{'trlen'} ? $info->{'trlen'} : "")."\t";
			print "".($info->{'cdslen'} ? $info->{'cdslen'} : "")."\t";
			print "".($info->{'exons'} ? $info->{'exons'} : "")."\t";
			print "".($info->{'description'} ? $info->{'description'} : "")."\t";
			
			print scalar(values(%{$genes{$p}{$s}{$g}{'all'}})), "\t";
			print scalar(values(%{$genes{$p}{$s}{$g}{'nonsyn'}})), "\t";
			my $first = 1;
			foreach my $v (keys(%{$genes{$p}{$s}{$g}{'all'}}))
			{
				print ";" if (!$first);
				print $v,"[",$genes{$p}{$s}{$g}{'all'}{$v},"]";
				$first = 0;
			}
			print "\n";
			
			$written ++;
		}		
	}
}
INFO("$written output lines written.");
