use warnings FATAL => qw( all );
use strict;

use List::Util qw(min max);

# STDIN: list with genes mutated in each patient (impacted-genes-list.tsv)
# STDOUT: gene/patient matrix indicating the number of times a gene is mutated across patients

use Getopt::Long;

# parse detailed results first
my ($mut_count, $mut_max_freq, $mut_details, $patient_ids);
GetOptions
(
	"mut-count" => \$mut_count,
	"mut-max-freq" => \$mut_max_freq,   
	"mut-details" => \$mut_details,   
	"patient-ids=s" => \$patient_ids   
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

# TABLE: impacted-genes.reseq.nonrel
my (%case_freq, %case_freq_ns, %case_freq_ns_af20, %mut_total, %mut_total_ns, %mut_gene_patient, %patients, %variants, %gene_info, %max_afs, %imp_exons, %imp_domains);

map { $patients{$_} = 1 } (split(",", $patient_ids)) if ($patient_ids);

<>; # skip header: patient\tcomparison\tgene\tchr\tstart\tend\ttr_len\tcds_len\texons\tcosmic\tdesc\tnum_mut\tnum_mut_nonsyn\tmax_af\tmax_af_ns\timp_exons\timp_exons_ns\tmut_effects\n
while(<>)
{
	chomp;
	my ($patient, $comp, $gene, $chr, $start, $end, $tr_len, $cds_len, $exons, $cosmic, $desc, $num_mutations, $num_mutations_nonsyn, $num_mutations_deleterious, $max_af, $max_af_ns, $ex, $ex_ns, $mutations, $domains) = split/\t/;
	
	$gene_info{$gene}{'chr'} = $chr;
	$gene_info{$gene}{'start'} = $start;
	$gene_info{$gene}{'end'} = $end;
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

	map { $imp_domains{$comp}{$gene}{$_} = 1 } split('\|', $domains); # impacted domains
		
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
		
		if ($max_af_ns >= 0.2)
		{
			if ($case_freq_ns_af20{$comp}{$gene})
			{
				$case_freq_ns_af20{$comp}{$gene} ++;
			}
			else
			{
				$case_freq_ns_af20{$comp}{$gene} = 1;
			}			
		}		
	}
}

my %all_genes;
foreach my $g (keys(%{$case_freq{'nonrel-dia'}}))
{
	$all_genes{$g} = 1;
}

my @sorted = sort { ($case_freq{'nonrel-dia'}{$b} ? $case_freq{'nonrel-dia'}{$b} : 0) <=> ($case_freq{'nonrel-dia'}{$a} ? $case_freq{'nonrel-dia'}{$a} : 0) } keys(%all_genes);

# TABLE: gene-patient-matrix
print "gene\tdescr\tchr\tstart\tend\texons\ttr_len\tcds_len\tcosmic\t";
print "freq\ttot\tfreq-ns\tfreq-ns-af20\ttot-ns\tmax-af\tmax-af-ns\timp-ex\timp-ex-ns\t";
map { print "$_\t" } keys(%patients);
print "\timp-domains";
print "\n";

foreach my $g (@sorted)
{
	print "$g\t";
	print $gene_info{$g}{'desc'},"\t",$gene_info{$g}{'chr'},"\t",$gene_info{$g}{'start'},"\t",$gene_info{$g}{'end'},"\t",$gene_info{$g}{'exons'},"\t",$gene_info{$g}{'tr_len'},"\t",$gene_info{$g}{'cds_len'},"\t",$gene_info{$g}{'cosmic'},"\t";

	print $case_freq{'nonrel-dia'}{$g} ? $case_freq{'nonrel-dia'}{$g} : "0", "\t";
	print $mut_total{'nonrel-dia'}{$g} ? $mut_total{'nonrel-dia'}{$g} : "0", "\t";
	print $case_freq_ns{'nonrel-dia'}{$g} ? $case_freq_ns{'nonrel-dia'}{$g} : "0", "\t";
	print $case_freq_ns_af20{'nonrel-dia'}{$g} ? $case_freq_ns_af20{'nonrel-dia'}{$g} : "0", "\t";
	print $mut_total_ns{'nonrel-dia'}{$g} ? $mut_total_ns{'nonrel-dia'}{$g} : "0", "\t";

	print $max_afs{'nonrel-dia'}{$g}{'all'} ? sprintf("%d", $max_afs{'nonrel-dia'}{$g}{'all'} * 100) : "", "\t";
	print $max_afs{'nonrel-dia'}{$g}{'ns'} ? sprintf("%d", $max_afs{'nonrel-dia'}{$g}{'ns'} * 100) : "", "\t";

	print join(",", sort { $a <=> $b } keys(%{$imp_exons{'nonrel-dia'}{$g}{'all'}})), "\t";
	print join(",", sort { $a <=> $b } keys(%{$imp_exons{'nonrel-dia'}{$g}{'ns'}})), "\t";

	foreach my $p (keys(%patients))
	{
		if ($mut_count)
		{
			print keys(%{$variants{'nonrel-dia'}{$g}{$p}}) > 0 ? scalar(keys(%{$variants{'nonrel-dia'}{$g}{$p}})) : " ", "\t";
		}
		elsif ($mut_max_freq)
		{
			print values(%{$variants{'nonrel-dia'}{$g}{$p}}) > 0 ? max(values(%{$variants{'nonrel-dia'}{$g}{$p}})) : " ", "\t";
		}
		else
		{
			print values(%{$variants{'nonrel-dia'}{$g}{$p}}) > 0 ? join("\|", values(%{$variants{'nonrel-dia'}{$g}{$p}})) : " ", "\t";
		}
	}

	print "\t", $imp_domains{'nonrel-dia'}{$g} ? join("|", keys(%{$imp_domains{'nonrel-dia'}{$g}})) : "";
	print "\n";
}
