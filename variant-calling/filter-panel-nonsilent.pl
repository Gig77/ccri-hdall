use warnings FATAL => qw( all );
use strict;

use Carp;

# read panel genes
my %pg;
open(PANEL, "/mnt/projects/hdall/results/panel-genes.tsv") or die "could not read panel gene file\n";
while(<PANEL>)
{
	chomp;
	$pg{$_} = 1;
}

while (<>)
{
	if (/^#/)
	{
		print $_;
		next;
	}
	
	chomp;
	my $line = $_;
	
	my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = split("\t", $line);
	
	next if ($id =~ /rs\d/); # benign common polymorphism
	next if ($filter eq "REM"); # present in remission
	
	my ($snpeff) = $info =~ /EFF=([^;\t]+)/;
	next if (!$snpeff);
		
	my ($gene, $effect) = get_impact($snpeff);
	next if (!$gene);
	
	my ($gt_t, $gq_t, $sdp_t, $dp_t, $rd_t, $ad_t, $freq_t, $pval_t, $rbq_t, $abq_t, $rdf_t, $rdr_t, $adf_t, $adr_t) = split(":", $sample);
	$freq_t =~ s/\%//; $freq_t = sprintf("%.2f", $freq_t/100);	
	
	print STDERR "$chrom\t$pos\t$ref\t$alt\t$freq_t\t$gene\t$effect\t$snpeff\n";
	print $line, "\n"; 	
}

sub get_impact
{
	my $effs = shift or die "ERROR: effect not specified";

	foreach my $eff (split(",", $effs))
	{
		my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
			or croak "ERROR: could not parse SNP effect: $effs";

		my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
			$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
				or croak "ERROR: could not parse SNP effect: $eff"; 

		next if (!$gene_name); # not affecting gene
		next if (!$pg{$gene_name}); # not affecting panel gene
		next if ($effect =~ /^(UPSTREAM|DOWNSTREAM|UTR_5_PRIME|UTR_5_DELETED|START_GAINED|UTR_3_PRIME|UTR_3_DELETED|INTRON|INTRON_CONSERVED|INTERGENIC|INTERGENIC_CONSERVED|INTRAGENIC|EXON|SYNONYMOUS_CODING)$/); # silent
		
		return ($gene_name, $effect);
	}
}
