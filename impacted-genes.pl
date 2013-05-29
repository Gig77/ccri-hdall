use warnings;
use strict;

use Data::Dumper;

my $debug = 1;
$| = 1; # turn on autoflush

# manually remap a few gene identifiers that do not match between SnpEff and UCSC, 
# presumably because minor version differences
my %remap = (
	'MST1P9' => 'MST1L',
	'FAM48B1' =>'SUPT20HL1',
	'HEATR7B2' => 'MROH2B',
	'ODZ4' => 'TENM4',
	'ODZ1' => 'TENM1',
	'USP17' => 'USP17L15',
	'CCDC165' => 'SOGA2',
	'C3orf15' => 'MAATS1',
	'C10orf140' => 'SKIDA1',
	'C19orf51' => 'DNAAF3',
	'PRIC285' => 'HELZ2',
	'FAM125B' => 'MVB12B',
	'BOD1L' => 'BOD1L1',
	'KDM4DL' => 'KDM4E',
	'FAM123B' => 'AMER1',
	'ZNF815' => 'ZNF815P',
	'LOC100192378' => 'ZFHX4-AS1',
	'HNRNPA1' => 'HNRNPA1P10',
	'CCDC13' => 'CCDC13-AS1',
	'MAGEA12' => 'CSAG4',
	'TXNDC5' => 'BLOC1S5-TXNDC5',
	'PGM5' => 'PGM5-AS1',
	'FLJ45983' => 'GATA3-AS1',
	'LOC643955' => 'ZNF733P',
	'LOC650623' => 'BEND3P3',
	
	'LOC84856' => 'LINC00839',
	'AK127292' => 'ANKRD20A19P',
	'LOC149837' => 'LINC00654',
	'LOC100132526' => 'FGD5P1',
	'LOC338739' => 'CSTF3-AS1',
	'TRAPPC2' => 'TRAPPC2P1',
	'AK302238' => 'GOLGA6L7P',
	'LOC154822' => 'LINC00689',
	'LOC126536' => 'LINC00661',
	'LOC646851' => 'FAM227A',
	'LOC649305' => 'LOC729732',
	'LOC146336' => 'SSTR5-AS1'

	# yet unmapped	
	# 'AHCTF1P1'
	# 'USP17L1P'
	# 'RNA45S5'
	# 'PTGER4P2-CDK2AP2P2'
);

## manually remap a few gene identifiers that do not match between SnpEff and UCSC ROI, 
## presumably because minor version differences
#my %remap = (
#	'MST1P9' => 'MST1L','MST1P9' => 'MST1L',
#	'FAM48B1' =>'SUPT20HL1','FAM48B1' =>'SUPT20HL1',
#	'HEATR7B2' => 'MROH2B','HEATR7B2' => 'MROH2B',
#	'ODZ4' => 'TENM4','ODZ4' => 'TENM4',
#	'ODZ1' => 'TENM1','ODZ1' => 'TENM1',
#	'USP17' => 'USP17L15','USP17' => 'USP17L15',
#	'CCDC165' => 'SOGA2','CCDC165' => 'SOGA2',
#	'C3orf15' => 'MAATS1','C3orf15' => 'MAATS1',
#	'C10orf140' => 'SKIDA1','C10orf140' => 'SKIDA1',
#	'C19orf51' => 'DNAAF3','C19orf51' => 'DNAAF3',
#	'PRIC285' => 'HELZ2','PRIC285' => 'HELZ2',
#	'FAM125B' => 'MVB12B','FAM125B' => 'MVB12B',
#	'BOD1L' => 'BOD1L1','BOD1L' => 'BOD1L1',
#	'KDM4DL' => 'KDM4E','KDM4DL' => 'KDM4E',
#	'FAM123B' => 'AMER1','FAM123B' => 'AMER1',
#	'ZNF815' => 'ZNF815P','ZNF815' => 'ZNF815P',
#	'LOC100192378' => 'ZFHX4-AS1','LOC100192378' => 'ZFHX4-AS1',
#	'SNRK' => 'SNRK-AS1','SNRK' => 'SNRK-AS1',
#	'HNRNPA1' => 'HNRNPA1P10',
#	'CCDC13' => 'CCDC13-AS1',
#	'MAGEA12' => 'CSAG4',
#	'TXNDC5' => 'BLOC1S5-TXNDC5',
#	'PGM5' => 'PGM5-AS1',
#	'FLJ45983' => 'GATA3-AS1',
#	'LOC643955' => 'ZNF733P',
#	'LOC650623' => 'BEND3P3'
#);

# read biomart id mapping
my (%prev2sym, %approved_symbols);
open(G, "$ENV{HOME}/hdall/results/gene-id-mapping.biomart-0.7.tsv") or die "ERROR: could not read gene list\n";
while(<G>)
{
	chomp;
	my ($approved_symbol, $entrez_gene_id, $accession_numbers, $approved_name, 
		$previous_symbols, $previous_names, $aliases, $name_aliases, $entrez_gene_id_ncbi) = split("\t");

	$approved_symbols{$approved_symbol} = 1;

	delete $prev2sym{$approved_symbol} # do not re-map approved symbols!
		if ($prev2sym{$approved_symbol}); 

	map {$prev2sym{$_} = $approved_symbol if (!$approved_symbols{$_})} split(", ", $previous_symbols);
	map {$prev2sym{$_} = $approved_symbol if (!$approved_symbols{$_})} split(", ", $aliases);
	map {$prev2sym{$_} = $approved_symbol if (!$approved_symbols{$_})} split(", ", $name_aliases);	
}
close(G);

# read gene information

my %id2sym;
my %sym2info;
open(G,"$ENV{HOME}/hdall/data/hg19/hg19.kgXref.txt") or die "could not open file $ENV{HOME}/hdall/data/hg19/hg19.kgXref.txt";
while(<G>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refSeq, $protAcc, $description, $rfamAcc, $tRnaName) = split(/\t/);

	next if (!$geneSymbol);
	$geneSymbol = $remap{$geneSymbol} if ($remap{$geneSymbol});	
	$geneSymbol = $prev2sym{$geneSymbol} if ($prev2sym{$geneSymbol});

	$id2sym{$kgID} = $geneSymbol;
	
	next if (!$description);

	#$sym2info{$prev2sym{$geneSymbol}}{'description'} = $description if ($prev2sym{$geneSymbol});
	$sym2info{$geneSymbol}{'description'} = $description;	
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
	#$geneSymbol = $remap{$geneSymbol} if ($remap{$geneSymbol});	
	#$geneSymbol = $prev2sym{$geneSymbol} if ($prev2sym{$geneSymbol});
	
	my $size = $chromEnd-$chromStart;
	
	# if multiple canonical transcripts for this gene symbol, use larger one
	next if ($canonical{$geneSymbol} and $sym2size{$geneSymbol} > $size); 
	
	$canonical{$geneSymbol} = $transcript;
	$sym2size{$geneSymbol} = $size;
}
close(G);

print STDERR "Canonical LDB3: ".$canonical{'LDB3'}."\n";

# manually add a few missing canonical transcripts
#$canonical{'uc002yit.1'} = 1;

open(G,"$ENV{HOME}/hdall/data/hg19/hg19.knownGene.txt") or die "could not open file $ENV{HOME}/hdall/data/hg19/hg19.knownGene.txt";
while(<G>)
{
	chomp;
	my ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $proteinID, $alignID) = split(/\t/);

	my $geneSymbol = $id2sym{$name};
	#$geneSymbol = $remap{$geneSymbol} if ($remap{$geneSymbol});	
	#$geneSymbol = $prev2sym{$geneSymbol} if ($prev2sym{$geneSymbol});
	
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
			if ($info->{'description'})
			{
				print "$info->{'description'}\t";
			}
			else
			{
				print STDERR "$0: WARNING: Could not map gene $g\n";
				print "\t";
			}
			
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
	
