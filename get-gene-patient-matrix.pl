use warnings FATAL => qw( all );
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

my (%case_freq, %case_freq_ns, %mut_total, %mut_total_ns, %mut_gene_patient, %patients, %variants, %gene_info, %max_afs, %imp_exons);
<>; # skip header
while(<>)
{
	chomp;
	my ($patient, $comp, $gene, $tr_len, $cds_len, $exons, $cosmic, $desc, $num_mutations, $num_mutations_nonsyn, $max_af, $max_af_ns, $ex, $ex_ns, $mutations) = split/\t/;
	
	$gene_info{$gene}{'tr_len'} = $tr_len;
	$gene_info{$gene}{'cds_len'} = $cds_len;
	$gene_info{$gene}{'exons'} = $exons;
	$gene_info{$gene}{'desc'} = $desc;
	$gene_info{$gene}{'cosmic'} = $cosmic;

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
		
		$max_af_ns = 0 if ($max_af_ns eq "");
		$max_afs{$comp}{$gene}{'all'} = $max_afs{$comp}{$gene}{'all'}
			? $max_afs{$comp}{$gene}{'all'} < $max_af ? $max_af : $max_afs{$comp}{$gene}{'all'}
			: $max_af;
		$max_afs{$comp}{$gene}{'ns'} = $max_afs{$comp}{$gene}{'ns'}
			? $max_afs{$comp}{$gene}{'ns'} < $max_af_ns ? $max_af_ns : $max_afs{$comp}{$gene}{'ns'}
			: $max_af_ns;
			
	} split(";", $mutations);
	
	map { $imp_exons{$comp}{$gene}{'all'}{$_} = 1 } split(",", $ex);  # impacted exons, all variants
	map { $imp_exons{$comp}{$gene}{'ns'}{$_} = 1 } split(",", $ex_ns); # impacted exons, only nonsynonymous variants
	
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

	if ($num_mutations_nonsyn > 0)
	{
		if ($case_freq_ns{$comp}{$gene})
		{
			$case_freq_ns{$comp}{$gene} ++;
		}
		else
		{
			$case_freq_ns{$comp}{$gene} = 1;
		}
	
		if ($mut_total_ns{$comp}{$gene})
		{
			$mut_total_ns{$comp}{$gene} += $num_mutations_nonsyn;
		}
		else
		{
			$mut_total_ns{$comp}{$gene} = $num_mutations_nonsyn;
		}		
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

print "gene\tdescr\texons\ttr_len\tcds_len\tcosmic\t";
print "freq-dia\ttot-dia\tfreq-dia-ns\ttot-dia-ns\tmax-af-dia\tmax-af-dia-ns\timp-ex-dia\timp-ex-dia-ns\t";
map { print "$_-dia\t" } keys(%patients);

print "freq-rel\ttot-rel\tfreq-rel-ns\ttot-rel-ns\tmax-af-rel\tmax-af-rel-ns\timp-ex-rel\timp-ex-rel-ns\t";
map { print "$_-rel\t" } keys(%patients);

print "freq-cons\t";
print "tot-cons";
map { print "\t$_-cons" } keys(%patients);
print "\n";

foreach my $g (@sorted)
{
	print "$g\t";
	print $gene_info{$g}{'desc'},"\t",$gene_info{$g}{'exons'},"\t",$gene_info{$g}{'tr_len'},"\t",$gene_info{$g}{'cds_len'},"\t",$gene_info{$g}{'cosmic'},"\t";

	print $case_freq{'rem_dia'}{$g} ? $case_freq{'rem_dia'}{$g} : "0", "\t";
	print $mut_total{'rem_dia'}{$g} ? $mut_total{'rem_dia'}{$g} : "0", "\t";
	print $case_freq_ns{'rem_dia'}{$g} ? $case_freq_ns{'rem_dia'}{$g} : "0", "\t";
	print $mut_total_ns{'rem_dia'}{$g} ? $mut_total_ns{'rem_dia'}{$g} : "0", "\t";

	print $max_afs{'rem_dia'}{$g}{'all'} ? sprintf("%d", $max_afs{'rem_dia'}{$g}{'all'} * 100) : "", "\t";
	print $max_afs{'rem_dia'}{$g}{'ns'} ? sprintf("%d", $max_afs{'rem_dia'}{$g}{'ns'} * 100) : "", "\t";

	print join(",", sort { $a <=> $b } keys(%{$imp_exons{'rem_dia'}{$g}{'all'}})), "\t";
	print join(",", sort { $a <=> $b } keys(%{$imp_exons{'rem_dia'}{$g}{'ns'}})), "\t";

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
	print $case_freq_ns{'rem_rel'}{$g} ? $case_freq_ns{'rem_rel'}{$g} : "0", "\t";
	print $mut_total_ns{'rem_rel'}{$g} ? $mut_total_ns{'rem_rel'}{$g} : "0", "\t";

	print $max_afs{'rem_rel'}{$g}{'all'} ? sprintf("%d", $max_afs{'rem_rel'}{$g}{'all'} * 100) : "", "\t";
	print $max_afs{'rem_rel'}{$g}{'ns'} ? sprintf("%d", $max_afs{'rem_rel'}{$g}{'ns'} * 100) : "", "\t";

	print join(",", sort { $a <=> $b } keys(%{$imp_exons{'rem_rel'}{$g}{'all'}})), "\t";
	print join(",", sort { $a <=> $b } keys(%{$imp_exons{'rem_rel'}{$g}{'ns'}})), "\t";

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
	print $mut_total{'cons'}{$g} ? $mut_total{'cons'}{$g} : "0";
	foreach my $p (keys(%patients))
	{
		print "\t";
		if ($mut_count)
		{
			print keys(%{$variants{'cons'}{$g}{$p}}) > 0 ? scalar(keys(%{$variants{'cons'}{$g}{$p}})) : " ";
		}
		elsif ($mut_max_freq)
		{
			print values(%{$variants{'cons'}{$g}{$p}}) > 0 ? max(values(%{$variants{'cons'}{$g}{$p}})) : " ";
		}
		else
		{
			print values(%{$variants{'cons'}{$g}{$p}}) > 0 ? join("\|", values(%{$variants{'cons'}{$g}{$p}})) : " ";
		}
	}
	print "\n";
}
