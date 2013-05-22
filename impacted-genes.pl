use warnings;
use strict;

my $debug = 1;
$| = 1; # turn on autoflush

my %genes;
<>; # skip header
while(<>)
{
	chomp;
	my ($patient, $sample, $var_type, $chr, $pos, $dbSNP, $ref, $alt, $gene, $add_genes, $impact, $effect, $depth_rem, $depth_leu, $freq, $snpeff) = split("\t");

	die "ERROR: $0: snpeff annotation missing from following line:\n$_\n"
		if (!$snpeff);

	next if ($effect =~ /(DOWNSTREAM|INTERGENIC|INTRON|UPSTREAM)/);
	
	$genes{$patient}{$sample}{$gene}{"$chr:$pos:$ref->$alt"} = $genes{$patient}{$sample}{$gene}{"$chr:$pos:$ref->$alt"} 
		? $genes{$patient}{$sample}{$gene}{"$chr:$pos:$ref->$alt"}.",$snpeff"
		: $snpeff;
}

print "patient\tcomparison\tgene\tnum_mut\tmut_effects\n";
foreach my $p (keys(%genes))
{
	foreach my $s (keys(%{$genes{$p}}))
	{
		my @sorted_genes = sort {values(%{$genes{$b}}) <=> values(%{$genes{$a}})} keys(%{$genes{$p}{$s}});
		foreach my $g (@sorted_genes)
		{
			print STDOUT $p, "\t", $s, "\t", $g, "\t", scalar(values(%{$genes{$p}{$s}{$g}})), "\t";
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
	
