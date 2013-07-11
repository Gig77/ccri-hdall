use warnings FATAL => qw( all );
use strict;

use Carp;
use Getopt::Long;

# STDIN: list with enriched pathways from DAVID
# STDOUT: pathway/patient matrix indicating the number of times a pathway is mutated across patients

my ($enriched_pathways_dia, $enriched_pathways_rel, $gene_patient_matrix);
GetOptions
(
	"enriched-pathways-dia=s" => \$enriched_pathways_dia,  
	"enriched-pathways-rel=s" => \$enriched_pathways_rel,  
	"gene-patient-matrix=s" => \$gene_patient_matrix  
);

my %enriched_pathways;

open(DIA,"$enriched_pathways_dia") or croak "ERROR: could not open file $enriched_pathways_dia\n"; 
<DIA>; # skip header
while(<DIA>)
{
	chomp;
	my ($pathway, $name, $class, $samples, $tot_var, $p_value, $fdr, $num_genes, $genes) = split /\t/;

	croak "ERROR: Could not parse pathway entry: $_\n" if (!$pathway or !$name or !$class or !$genes);
	$enriched_pathways{"$pathway\t$name\t$class"}{dia} = "$p_value\t$genes";
}
close(DIA);

open(REL,"$enriched_pathways_rel") or croak "ERROR: could not open file $enriched_pathways_rel\n"; 
<REL>; # skip header
while(<REL>)
{
	chomp;
	my ($pathway, $name, $class, $samples, $tot_var, $p_value, $fdr, $num_genes, $genes) = split /\t/;

	croak "ERROR: Could not parse pathway entry: $_\n" if (!$pathway or !$name or !$class or !$genes);
	$enriched_pathways{"$pathway\t$name\t$class"}{rel} = "$p_value\t$genes";
}
close(REL);

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
print "Pathway\tName\tClass\tp-dia\tfreq-dia";
map { print "\t$_" } (@patients_dia);
print "\tp-rel\tfreq-rel";
map { print "\t$_" } (@patients_rel);
print "\tGenes.dia\tGenes.rel\n";

foreach my $p (keys(%enriched_pathways))
{
	my ($pathway, $name, $class) = split("\t", $p);

	print "$pathway\t$name\t$class\t";

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