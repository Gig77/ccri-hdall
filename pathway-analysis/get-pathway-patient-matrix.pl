use warnings FATAL => qw( all );
use strict;

use Carp;
use Getopt::Long;
use List::Util qw(min max);

# STDIN: list with enriched pathways from DAVID
# STDOUT: pathway/patient matrix indicating the number of times a pathway is mutated across patients

my ($enriched_pathways_dia, $enriched_pathways_rel, $curated_pathways_file, $filtered_variants);
GetOptions
(
	"enriched-pathways-dia=s" => \$enriched_pathways_dia,  
	"enriched-pathways-rel=s" => \$enriched_pathways_rel,  
	"curated-pathways=s" => \$curated_pathways_file, # manually curated pathways with recurrently mutated genes at relapse 
	"filtered-variants=s" => \$filtered_variants  
);

# TABLE: filtered-variants
my %variants;
open(VAR,"$filtered_variants") or croak "ERROR: could not open file $filtered_variants\n"; 
<VAR>; # skip header
while(<VAR>)
{
	chomp;
	my ($patient, $sample, $var_type, $status, $chr, $pos, $dbSNP, $ref, $alt, $gene, $add_genes, $impact, $effect, $exons, 
		$dp_rem_tot, $dp_rem_ref, $dp_rem_var, $freq_rem, $dp_leu_tot, $dp_leu_ref, $dp_leu_var, $freq_leu) = split("\t");

	next if ($effect =~ /^(DOWNSTREAM|INTERGENIC|INTRON|UPSTREAM|INTERGENIC_CONSERVED|SYNONYMOUS_START|SYNONYMOUS_CODING|SYNONYMOUS_STOP|UTR_5_PRIME|UTR_5_DELETED|START_GAINED|UTR_3_PRIME|UTR_3_DELETED|INTRON_CONSERVED|INTRAGENIC|EXON)$/);
	
	$variants{$patient}{$sample}{$gene} = $variants{$patient}{$sample}{$gene} ? max($freq_leu, $variants{$patient}{$sample}{$gene}) : $freq_leu; 
}
close(VAR);

# TABLE: sm_pathways.annotated
my %enriched_pathways;
open(DIA,"$enriched_pathways_dia") or croak "ERROR: could not open file $enriched_pathways_dia\n"; 
<DIA>; # skip header
while(<DIA>)
{
	chomp;
	my ($pathway, $name, $class, $size, $samples, $tot_var, $p_value, $fdr, $num_genes, $genes, $link) = split /\t/;

	croak "ERROR: Could not parse pathway entry: $_\n" if (!$pathway or !$name or !$class or !$genes);
	$enriched_pathways{"$pathway\t$name\t$class\t$size"}{dia} = "$p_value\t$genes";
}
close(DIA);

# TABLE: sm_pathways.annotated
open(REL,"$enriched_pathways_rel") or croak "ERROR: could not open file $enriched_pathways_rel\n"; 
<REL>; # skip header
while(<REL>)
{
	chomp;
	my ($pathway, $name, $class, $size, $samples, $tot_var, $p_value, $fdr, $num_genes, $genes, $link) = split /\t/;

	croak "ERROR: Could not parse pathway entry: $_\n" if (!$pathway or !$name or !$class or !$genes);
	$enriched_pathways{"$pathway\t$name\t$class\t$size"}{rel} = "$p_value\t$genes";
}
close(REL);

my %curated_pathways;
open(CUR,"$curated_pathways_file") or croak "ERROR: could not open file $curated_pathways_file\n"; 
<CUR>; # skip header
while(<CUR>)
{
	chomp;
	my ($gene, $group1, $group2, $group3, $category, $comment, $source) = split /\t/;
	die "$_\n" if ($group1 eq "Renate");

	croak "ERROR: coud not parse following line:\n$_\n" if (!$group1 or !$gene);
	$category = "" if (!$category);
	
	if ($curated_pathways{"$group1\t$group1\t$category"})
	{
		$curated_pathways{"$group1\t$group1\t$category"} = $curated_pathways{"$group1\t$group1\t$category"}.",".$gene;
	}
	else
	{
		$curated_pathways{"$group1\t$group1\t$category"} = "$gene";	
	}
	
	if ($group2)
	{
		if ($curated_pathways{"$group2\t$group2\t$category"})
		{
			$curated_pathways{"$group2\t$group2\t$category"} = $curated_pathways{"$group2\t$group2\t$category"}.",".$gene;
		}
		else
		{
			$curated_pathways{"$group2\t$group2\t$category"} = "$gene";	
		}
	}

	if ($group3)
	{
		if ($curated_pathways{"$group3\t$group3\t$category"})
		{
			$curated_pathways{"$group3\t$group3\t$category"} = $curated_pathways{"$group3\t$group3\t$category"}.",".$gene;
		}
		else
		{
			$curated_pathways{"$group3\t$group3\t$category"} = "$gene";	
		}
	}
}
close(CUR);

# output header
# TABLE: pathway-patient-matrix
print "Pathway\tName\tClass\tSize\tp-dia\tfreq-dia";
map { print "\t$_-dia" } (keys(%variants));
print "\tp-rel\tfreq-rel";
map { print "\t$_-rel" } (keys(%variants));
print "\tGenes.dia\tGenes.rel\n";

# music enriched pathways
foreach my $p (keys(%enriched_pathways))
{
	my ($pathway, $name, $class, $size) = split("\t", $p);

	print "$pathway\t$name\t$class\t$size\t";

	my ($genes_dia, $genes_rel) = ("", "");
	if ($enriched_pathways{$p}{dia})
	{
		my ($p_value, $genes) = split("\t", $enriched_pathways{$p}{dia});
		$genes_dia = $genes;
		my ($line, $tot_num_mut) = ("", 0);
		foreach my $pd (keys(%variants))
		{
			my @mut_genes;
			foreach my $g (split(",", $genes))
			{
				$g =~ s/\(\d+\)//;
				if ($variants{$pd}{'rem_dia'}{$g})
				{
					push(@mut_genes, "$g(".$variants{$pd}{'rem_dia'}{$g}.")");
				}
			}
			$line .= "\t".(@mut_genes > 0 ? join(",", @mut_genes) : " ");
			$tot_num_mut ++ if (@mut_genes > 0);
		}
		print "$p_value\t$tot_num_mut$line"
	}
	else
	{
		print "1\t0";
		map { print "\t" } (keys(%variants));
	}
	print "\t";		
	
	if ($enriched_pathways{$p}{rel})
	{
		my ($p_value, $genes) = split("\t", $enriched_pathways{$p}{rel});
		$genes_rel = $genes;
		my ($line, $tot_num_mut) = ("", 0);
		foreach my $pr (keys(%variants))
		{
			my @mut_genes;
			foreach my $g (split(",", $genes))
			{
				$g =~ s/\(\d+\)//;
				if ($variants{$pr}{'rem_rel'}{$g})
				{
					push(@mut_genes, "$g(".$variants{$pr}{'rem_rel'}{$g}.")");
				}
			}
			$line .= "\t".(@mut_genes > 0 ? join(",", @mut_genes) : " ");
			$tot_num_mut ++ if (@mut_genes > 0);
		}
		print "$p_value\t$tot_num_mut$line"		
	}
	else
	{
		print "1\t0";
		map { print "\t" } (keys(%variants));
	}
	print "\t$genes_dia\t$genes_rel\n";
}

# curated pathways with recurrently mutated genes
foreach my $p (keys(%curated_pathways))
{
	my ($pathway, $name, $class) = split("\t", $p);
	my $genes = $curated_pathways{$p};
	my $outputline = "";

	$outputline .= "$pathway".($class ? ":$class" : "")."\t";
	$outputline .= "$name".($class ? ":$class" : "")."\t";
	$outputline .= "CURATED\t";
	$outputline .= "\t"; # put size here

	my %genes_dia;
	{
		my ($line, $tot_num_mut) = ("", 0);
		foreach my $pd (keys(%variants))
		{
			my @mut_genes;
			foreach my $g (split(",", $genes))
			{
				$g =~ s/\(\d+\)//;
				if ($variants{$pd}{'rem_dia'}{$g})
				{
					push(@mut_genes, "$g(".$variants{$pd}{'rem_dia'}{$g}.")");
					$genes_dia{$g} = 1;
				}
			}
			$line .= "\t".(@mut_genes > 0 ? join(",", @mut_genes) : " ");
			$tot_num_mut ++ if (@mut_genes > 0);
		}
		$outputline .= "1\t$tot_num_mut$line"
	}
	$outputline .= "\t";		
	
	my %genes_rel;
	{
		my ($line, $tot_num_mut) = ("", 0);
		foreach my $pr (keys(%variants))
		{
			my @mut_genes;
			foreach my $g (split(",", $genes))
			{
				$g =~ s/\(\d+\)//;
				if ($variants{$pr}{'rem_rel'}{$g})
				{
					push(@mut_genes, "$g(".$variants{$pr}{'rem_rel'}{$g}.")");
					$genes_rel{$g} = 1;
				}
			}
			$line .= "\t".(@mut_genes > 0 ? join(",", @mut_genes) : " ");
			$tot_num_mut ++ if (@mut_genes > 0);
		}
		$outputline .= "1\t$tot_num_mut$line"		
	}

	$outputline .= "\t".join(",", keys(%genes_dia))."\t".join(",", keys(%genes_rel))."\n";
	
	print $outputline 
		if (keys(%genes_dia) > 0 or keys(%genes_rel) > 0);
}