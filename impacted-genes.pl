use warnings;
use strict;

my $debug = 1;
$| = 1; # turn on autoflush

my %genes;
<>; # skip header
while(<>)
{
	chomp;
	my ($patient, $sample, $var_type, $chr, $pos, $dbSNP, $ref, $alt, $gene, $impact, $variant_status, $depth_rem, $depth_leu, $freq, $effects) = split("\t");

	foreach my $eff (split(",", $effects))
	{
		my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
			or die "could not parse SNP effect: $chr:$pos:$eff\n";
		next if ($effect =~ /(DOWNSTREAM|INTERGENIC|INTRON|UPSTREAM)/);

		my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
			$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
				or die "could not parse SNP effect: $chr:$pos:$rest\n";

#		print "$patient\t$sample\t$gene_name\t$chr:$pos:$ref->$alt\n";
		$genes{$patient}{$sample}{$gene_name}{"$chr:$pos:$ref->$alt"} = $genes{$patient}{$sample}{$gene_name}{"$chr:$pos:$ref->$alt"} 
			? $genes{$patient}{$sample}{$gene_name}{"$chr:$pos:$ref->$alt"}.",$eff"
			: $eff;
	}
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
	
