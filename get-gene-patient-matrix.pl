use warnings;
use strict;

# STDIN: list with genes mutated in each patient (impacted-genes-list.tsv)
# STDOUT: gene/patient matrix indicating the number of times a gene is mutated across patients

my (%case_freq, %mut_total, %mut_gene_patient, %patients, %variants);
<>; # skip header
while(<>)
{
	chomp;
	my ($patient, $comp, $gene, $num_mutations, $mutations) = split/\t/;

	$patients{$patient} = 1;
	map { $variants{$comp}{$gene}{$patient}{$_} = 1 } split(";", $mutations);
	$mut_gene_patient{$comp}{$gene}{$patient} = $num_mutations;

	if ($case_freq{$comp}{$gene})
	{
		$case_freq{$comp}{$gene} ++;
	}
	else
	{
		$case_freq{$comp}{$gene} = 1;
	}

	if ($mut_total{$comp}{$gene})
	{
		$mut_total{$comp}{$gene} += $num_mutations;
	}
	else
	{
		$mut_total{$comp}{$gene} = $num_mutations;
	}
}

#my (%freq_cons, %mut_cons);
# find mutations conserved b/w diagnosis and relapse
my %all_genes;
foreach my $g (keys(%{$case_freq{'rem_dia'}}))
{
	$all_genes{$g} = 1;
	
	foreach my $p (keys(%patients))
	{
		my $counted = 0;
		foreach my $v (keys(%{$variants{'rem_dia'}{$g}{$p}}))
		{
			if (exists $variants{'rem_rel'}{$g}{$p}{$v})  # variant also in relapse?
			{
				$mut_gene_patient{'cons'}{$g}{$p} = $mut_gene_patient{'cons'}{$g}{$p} ? $mut_gene_patient{'cons'}{$g}{$p} + 1 : 1;  
				$mut_total{'cons'}{$g} = $mut_total{'cons'}{$g} ? $mut_total{'cons'}{$g} + 1 : 1;
				$case_freq{'cons'}{$g} = $case_freq{'cons'}{$g} ? $case_freq{'cons'}{$g} + 1 : 1
					if (!$counted);
				$counted = 1;
			}
#			print "$g\t$p\t$v\n";
		}
	}
}
foreach my $g (keys(%{$case_freq{'rem_rel'}}))
{
	$all_genes{$g} = 1;
}

my @sorted = sort { ($case_freq{'cons'}{$b} ? $case_freq{'cons'}{$b} : 0) <=> ($case_freq{'cons'}{$a} ? $case_freq{'cons'}{$a} : 0) } keys(%all_genes);

print "\tfreq-dia\t";
print "tot-dia\t";
map { print "$_-dia\t" } keys(%patients);

print "freq-rel\t";
print "tot-rel\t";
map { print "$_-rel\t" } keys(%patients);

print "freq-cons\t";
print "tot-cons";
map { print "\t$_-cons" } keys(%patients);
print "\n";

foreach my $g (@sorted)
{
	print "$g\t";

	print $case_freq{'rem_dia'}{$g} ? $case_freq{'rem_dia'}{$g} : "", "\t";
	print $mut_total{'rem_dia'}{$g} ? $mut_total{'rem_dia'}{$g} : "", "\t";
	foreach my $p (keys(%patients))
	{
		print $mut_gene_patient{'rem_dia'}{$g}{$p} ? $mut_gene_patient{'rem_dia'}{$g}{$p} : "", "\t";
	}

	print $case_freq{'rem_rel'}{$g} ? $case_freq{'rem_rel'}{$g} : "", "\t";
	print $mut_total{'rem_rel'}{$g} ? $mut_total{'rem_rel'}{$g} : "", "\t";
	foreach my $p (keys(%patients))
	{
		print $mut_gene_patient{'rem_rel'}{$g}{$p} ? $mut_gene_patient{'rem_rel'}{$g}{$p} : "", "\t";
	}

	print $case_freq{'cons'}{$g} ? $case_freq{'cons'}{$g} : "", "\t";
	print $mut_total{'cons'}{$g} ? $mut_total{'cons'}{$g} : "";
	foreach my $p (keys(%patients))
	{
		print "\t", $mut_gene_patient{'cons'}{$g}{$p} ? $mut_gene_patient{'cons'}{$g}{$p} : "";
	}
	print "\n";
}
