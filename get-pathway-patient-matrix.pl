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
	my ($category, $term, $count, $count_perc, $p_value, $genes, $list_total, $pop_hits, $pop_total, $fold_enrichment, $bonferroni, $benjamini, $fdr) = split /\t/;

	$enriched_pathways{"$category\t$term"}{dia} = "$p_value\t$genes";
}
close(DIA);

open(REL,"$enriched_pathways_rel") or croak "ERROR: could not open file $enriched_pathways_rel\n"; 
<REL>; # skip header
while(<REL>)
{
	chomp;
	my ($category, $term, $count, $count_perc, $p_value, $genes, $list_total, $pop_hits, $pop_total, $fold_enrichment, $bonferroni, $benjamini, $fdr) = split /\t/;

	$enriched_pathways{"$category\t$term"}{rel} = "$p_value\t$genes";
}
close(REL);

my %gene_patient;
open(MATRIX, "$gene_patient_matrix") or croak "ERROR: could not open file $gene_patient_matrix\n";

my (@patients_dia, @patients_rel);
my $header = <MATRIX>;
chomp($header);
my @hfields = split("\t", $header);
for (my $d = 3; $d <= 22; $d ++) { push(@patients_dia, $hfields[$d]); }
for (my $d = 25; $d <= 44; $d ++) { push(@patients_rel, $hfields[$d]); }

while(<MATRIX>)
{
	chomp;
	my @fields = split /\t/;
	my $gene = $fields[0];
	
	for (my $d = 3; $d <= 22; $d ++) { $gene_patient{$gene}{$patients_dia[$d-3]} = $fields[$d]; }
	for (my $d = 25; $d <= 44; $d ++) { $gene_patient{$gene}{$patients_rel[$d-25]} = $fields[$d]; }
}
close(MATRIX);

# output header
print "category\tterm\tp-dia\tfreq-dia";
map { print "\t$_" } (@patients_dia);
print "\tp-rel\tfreq-rel";
map { print "\t$_" } (@patients_rel);
print "\tgenes-dia\tgenes-rel\n";

foreach my $p (keys(%enriched_pathways))
{
	my ($category, $term) = split("\t", $p);

	print "$category\t$term\t";

	my ($genes_dia, $genes_rel) = ("", "");
	if ($enriched_pathways{$p}{dia})
	{
		my ($p_value, $genes) = split("\t", $enriched_pathways{$p}{dia});
		$genes_dia = $genes;
		my ($line, $tot_num_mut) = ("", 0);
		foreach my $pd (@patients_dia)
		{
			my $num_mut = 0;
			foreach my $g (split(",", $genes))
			{
				$num_mut ++ if ($gene_patient{$g}{$pd});
			}
			$line .= "\t".($num_mut > 0 ? $num_mut : "");
			$tot_num_mut ++ if ($num_mut > 0);
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
			my $num_mut = 0;
			foreach my $g (split(",", $genes))
			{
				$num_mut ++ if ($gene_patient{$g}{$pr});
			}
			$line .= "\t".($num_mut > 0 ? $num_mut : "");
			$tot_num_mut ++ if ($num_mut > 0);
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