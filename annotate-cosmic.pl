use warnings FATAL => qw( all );
use strict;

use lib "$ENV{HOME}/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Carp;
use Getopt::Long;

# parse detailed results first
my ($cosmic_mutation_file);
GetOptions
(
	"cosmic-mutation-file=s" => \$cosmic_mutation_file
);

croak "ERROR: --cosmic-mutation-file not specified" if (!$cosmic_mutation_file);

# read cosmic mutations
my %cosmic;
my $entries_read = 0;
open(C, "$cosmic_mutation_file") or croak "ERROR: Could not open file $cosmic_mutation_file";
<C>; # skip header
while(<C>)
{
	chomp;
	my ($gene_name, $accession_number, $gene_cds_length, $hgnc_id, $sample_name, $id_sample, $id_tumour, $primary_site, $site_subtype, $primary_histology,
		$histology_subtype, $genome_wide_screen, $mutation_id, $mutation_cds, $mutation_aa, $mutation_description, $mutation_zygosity, $mutation_ncbi36_genome_position,
		$mutation_ncbi36_strand, $mutation_GRCh37_genome_position, $mutation_GRCh37_strand, $mutation_somatic_status, $pubmed_pmid, $sample_source, 
		$tumour_origin, $age, $comments) = split /\t/;
	
	$cosmic{$mutation_GRCh37_genome_position} = defined $cosmic{$mutation_GRCh37_genome_position} ? $cosmic{$mutation_GRCh37_genome_position} + 1 : 1;
	
	if ($mutation_aa =~ /p\.(.)(\d+)(.+)/)
	{
		my ($prev_aa, $aa_pos, $after_aa) = ($1, $2, $3);
		$cosmic{"$gene_name:$prev_aa:$aa_pos"} = defined $cosmic{"$gene_name:$prev_aa:$aa_pos"} ? $cosmic{"$gene_name:$prev_aa:$aa_pos"} + 1 : 1;
	}
	$entries_read ++;
}
close(C);
INFO("$entries_read mutations read from file $cosmic_mutation_file");

# TABLE: filtered-variants.cosmic
# TABLE: filtered-variants
my $header = <>;
chomp($header);
print "$header\tcosmic_hits_nt\tcosmic_hits_aa\n";
while(<>)
{
	chomp;
	my ($patient, $sample, $var_type, $status, $chr, $pos, $dbSNP, $ref, $alt, $gene, $add_genes, $impact, $effect, $exons, 
		$dp_rem_tot, $dp_rem_ref, $dp_rem_var, $freq_rem, $dp_leu_tot, $dp_leu_ref, $dp_leu_var, $freq_leu, $aa_change, $snpeff,
		$polyphen2, $sift, $gerp, $siphy, $interpro, $af_1000g) = split("\t");

	my $loc = "$chr:$pos-$pos";
	$loc =~ s/^chr//;
	
	print "$_\t";
	print $cosmic{$loc} ? $cosmic{$loc} : "0", "\t";
	print aa_hits([$gene, split(",", $add_genes)], $snpeff), "\n";	
}

# -------

sub aa_hits
{
	my $genes = shift;
	my $snpeff = shift;
	
	return "non-coding" if (@$genes == 0);
	
	my $aa_change_found = 0; 
	foreach my $gene (@$genes) # check each gene
	{
		foreach my $eff (split(",", $snpeff)) # check all isoforms for cosmic match
		{
			my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
				or croak "ERROR: could not parse SNP effect: $snpeff";
	
			my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
				$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
					or croak "ERROR: could not parse SNP effect: $eff";
					 
			if ($aa_change =~ /(.)(\d+)(.+)/)
			{
				$aa_change_found = 1;			
				my ($prev_aa, $aa_pos, $after_aa) = ($1, $2, $3);
				return $cosmic{"$gene:$prev_aa:$aa_pos"} if (defined $cosmic{"$gene:$prev_aa:$aa_pos"});
			}
		}				
	}
	
	return "non-coding" if (!$aa_change_found);
	return "0";
}