use warnings FATAL => qw( all );
use strict;

use Carp;
use Getopt::Long;

# STDIN: list with enriched pathways from DAVID
# STDOUT: pathway/patient matrix indicating the number of times a pathway is mutated across patients

my ($enriched_pathways_dia, $enriched_pathways_rel, $curated_pathways_file, $gene_patient_matrix);
GetOptions
(
	"enriched-pathways-dia=s" => \$enriched_pathways_dia,  
	"enriched-pathways-rel=s" => \$enriched_pathways_rel,  
	"curated-pathways=s" => \$curated_pathways_file, # manually curated pathways with recurrently mutated genes at relapse 
	"gene-patient-matrix=s" => \$gene_patient_matrix  
);

my %enriched_pathways;

# TABLE: sm_pathways.annotated
open(DIA,"$enriched_pathways_dia") or croak "ERROR: could not open file $enriched_pathways_dia\n"; 
<DIA>; # skip header
while(<DIA>)
{
	chomp;
	my ($pathway, $name, $class, $size, $samples, $tot_var, $p_value, $fdr, $num_genes, $genes) = split /\t/;

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
	my ($pathway, $name, $class, $size, $samples, $tot_var, $p_value, $fdr, $num_genes, $genes) = split /\t/;

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

# TABLE: gene-patient-matrix
my %gene_patient;
open(MATRIX, "$gene_patient_matrix") or croak "ERROR: could not open file $gene_patient_matrix\n";

my (@patients_dia, @patients_rel);
my $header = <MATRIX>;
chomp($header);
my @hfields = split("\t", $header);
print STDERR "Patients diagnosis:";
for (my $d = 17; $d <= 36; $d ++) { print STDERR " $hfields[$d]"; push(@patients_dia, $hfields[$d]); }
print STDERR "\nPatients relapse:";
for (my $d = 45; $d <= 64; $d ++) { print STDERR " $hfields[$d]"; push(@patients_rel, $hfields[$d]); }
print STDERR "\n";

while(<MATRIX>)
{
	chomp;
	my @fields = split /\t/;
	my $gene = $fields[0];
	
	for (my $d = 17; $d <= 36; $d ++) { $gene_patient{$gene}{$patients_dia[$d-17]} = $fields[$d]; }
	for (my $d = 45; $d <= 64; $d ++) { $gene_patient{$gene}{$patients_rel[$d-45]} = $fields[$d]; }
}
close(MATRIX);

# output header
# TABLE: pathway-patient-matrix
print "Pathway\tName\tClass\tSize\tp-dia\tfreq-dia";
map { print "\t$_" } (@patients_dia);
print "\tp-rel\tfreq-rel";
map { print "\t$_" } (@patients_rel);
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
		foreach my $pd (@patients_dia)
		{
			my @mut_genes;
			foreach my $g (split(",", $genes))
			{
				if ($gene_patient{$g}{$pd} and $gene_patient{$g}{$pd} ne " ")
				{
					push(@mut_genes, $g);
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
		map { print "\t" } (@patients_dia);
	}
	print "\t";		
	
	if ($enriched_pathways{$p}{rel})
	{
		my ($p_value, $genes) = split("\t", $enriched_pathways{$p}{rel});
		$genes_rel = $genes;
		my ($line, $tot_num_mut) = ("", 0);
		foreach my $pr (@patients_rel)
		{
			my @mut_genes;
			foreach my $g (split(",", $genes))
			{
				if ($gene_patient{$g}{$pr} and $gene_patient{$g}{$pr} ne " ")
				{
					push(@mut_genes, $g);					
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
		map { print "\t" } (@patients_rel);
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
		foreach my $pd (@patients_dia)
		{
			my @mut_genes;
			foreach my $g (split(",", $genes))
			{
				if ($gene_patient{$g}{$pd} and $gene_patient{$g}{$pd} ne " ")
				{
					push(@mut_genes, $g);
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
		foreach my $pr (@patients_rel)
		{
			my @mut_genes;
			foreach my $g (split(",", $genes))
			{
				if ($gene_patient{$g}{$pr} and $gene_patient{$g}{$pr} ne " ")
				{
					push(@mut_genes, $g);					
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