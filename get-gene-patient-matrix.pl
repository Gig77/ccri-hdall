use warnings;
use strict;

use List::Util qw(min max);

# STDIN: list with genes mutated in each patient (impacted-genes-list.tsv)
# STDOUT: gene/patient matrix indicating the number of times a gene is mutated across patients

use Getopt::Long;

# parse detailed results first
my ($mut_count, $mut_max_freq, $mut_details);
GetOptions
(
	"mut-count" => \$mut_count,
	"mut-max-freq" => \$mut_max_freq,   
	"mut-details" => \$mut_details   
);

my %impact2flag =
(
	'HIGH' => 'H',
	'MODERATE' => 'M',
	'LOW' => 'L',
	'MODIFIER' => 'O'
);

die "ERROR: invalid or missing list type\n"
	if (!$mut_count and !$mut_max_freq and !$mut_details);

my (%case_freq, %mut_total, %mut_gene_patient, %patients, %variants, %gene_info);
<>; # skip header
while(<>)
{
	chomp;
	my ($patient, $comp, $gene, $tr_len, $cds_len, $exons, $desc, $num_mutations, $mutations) = split/\t/;
	
	$gene_info{$gene}{'tr_len'} = $tr_len;
	$gene_info{$gene}{'cds_len'} = $cds_len;
	$gene_info{$gene}{'exons'} = $exons;
	$gene_info{$gene}{'desc'} = $desc;

	$patients{$patient} = 1;
	map { 
		my ($chr, $start, $change, $freq, $impact, $effect) = split(":");
		
		if ($mut_details)
		{
			$variants{$comp}{$gene}{$patient}{"$chr:$start:$change"} = sprintf("%d(%s)", $freq*100, $impact2flag{$impact}); 	
		}
		else
		{
			$variants{$comp}{$gene}{$patient}{"$chr:$start:$change"} = sprintf("%d", $freq*100); 			
		}
	} split(";", $mutations);
	
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
				my ($freq_dia) = $variants{'rem_dia'}{$g}{$p}{$v} =~ /(\d+)/;
				my ($freq_rel) = $variants{'rem_rel'}{$g}{$p}{$v} =~ /(\d+)/;
				
				$variants{'cons'}{$g}{$p}{$v} = $mut_max_freq ? max($freq_dia, $freq_rel) : "$freq_dia>$freq_rel";  
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

print "gene\tdescr\texons\ttr_len\tcds_len\tfreq-dia\t";
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
	print $gene_info{$g}{'desc'},"\t",$gene_info{$g}{'exons'},"\t",$gene_info{$g}{'tr_len'},"\t",$gene_info{$g}{'cds_len'},"\t";

	print $case_freq{'rem_dia'}{$g} ? $case_freq{'rem_dia'}{$g} : "0", "\t";
	print $mut_total{'rem_dia'}{$g} ? $mut_total{'rem_dia'}{$g} : "0", "\t";
	foreach my $p (keys(%patients))
	{
		if ($mut_count)
		{
			print keys(%{$variants{'rem_dia'}{$g}{$p}}) > 0 ? scalar(keys(%{$variants{'rem_dia'}{$g}{$p}})) : " ", "\t";
		}
		elsif ($mut_max_freq)
		{
			print values(%{$variants{'rem_dia'}{$g}{$p}}) > 0 ? max(values(%{$variants{'rem_dia'}{$g}{$p}})) : " ", "\t";
		}
		else
		{
			print values(%{$variants{'rem_dia'}{$g}{$p}}) > 0 ? join("\|", values(%{$variants{'rem_dia'}{$g}{$p}})) : " ", "\t";
		}
	}

	print $case_freq{'rem_rel'}{$g} ? $case_freq{'rem_rel'}{$g} : "0", "\t";
	print $mut_total{'rem_rel'}{$g} ? $mut_total{'rem_rel'}{$g} : "0", "\t";
	foreach my $p (keys(%patients))
	{
		if ($mut_count)
		{
			print keys(%{$variants{'rem_rel'}{$g}{$p}}) > 0 ? scalar(keys(%{$variants{'rem_rel'}{$g}{$p}})) : " ", "\t";
		}
		elsif ($mut_max_freq)
		{
			print values(%{$variants{'rem_rel'}{$g}{$p}}) > 0 ? max(values(%{$variants{'rem_rel'}{$g}{$p}})) : " ", "\t";
		}
		else
		{
			print values(%{$variants{'rem_rel'}{$g}{$p}}) > 0 ? join("\|", values(%{$variants{'rem_rel'}{$g}{$p}})) : " ", "\t";
		}
	}

	print $case_freq{'cons'}{$g} ? $case_freq{'cons'}{$g} : "0", "\t";
	print $mut_total{'cons'}{$g} ? $mut_total{'cons'}{$g} : "0", "\t";
	foreach my $p (keys(%patients))
	{
		if ($mut_count)
		{
			print keys(%{$variants{'cons'}{$g}{$p}}) > 0 ? scalar(keys(%{$variants{'cons'}{$g}{$p}})) : " ", "\t";
		}
		elsif ($mut_max_freq)
		{
			print values(%{$variants{'cons'}{$g}{$p}}) > 0 ? max(values(%{$variants{'cons'}{$g}{$p}})) : " ", "\t";
		}
		else
		{
			print values(%{$variants{'cons'}{$g}{$p}}) > 0 ? join("\|", values(%{$variants{'cons'}{$g}{$p}})) : " ", "\t";
		}
	}
	print "\n";
}
