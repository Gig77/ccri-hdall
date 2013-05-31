use warnings;
use strict;

use lib "$ENV{HOME}/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Carp;
use Getopt::Long;

# parse detailed results first
my ($gene_patient_matrix, $smg_dia, $smg_rel, $smp_dia, $smp_rel);
GetOptions
(
	"gene-patient-matrix=s" => \$gene_patient_matrix,
	"smg-dia=s" => \$smg_dia,   
	"smg-rel=s" => \$smg_rel,   
	"smp-dia=s" => \$smp_dia,  
	"smp-rel=s" => \$smp_rel   
);

# read significantly mutated genes at diagnosis
my %smg_dia_pvalue;
open(D,"$smg_dia") or croak "ERROR: could not read file $smg_dia\n";
<D>; # skip header
while(<D>)
{
	chomp;
	my ($gene, $indels, $snvs, $totmut, $covd_bps, $mut_per_mb, $p_fcpt, $p_lrt, $p_ct, $fdr_fcpt, $fdr_lrt, $fdr_ct) = split(/\t/);
	$smg_dia_pvalue{$gene} = $p_ct;
}
close(D);
INFO(scalar(keys(%smg_dia_pvalue))." genes read from file $smg_dia");

# read significantly mutated genes at relapse
my %smg_rel_pvalue;
open(R,"$smg_rel") or croak "ERROR: could not read file $smg_rel\n";
<R>; # skip header
while(<R>)
{
	chomp;
	my ($gene, $indels, $snvs, $totmut, $covd_bps, $mut_per_mb, $p_fcpt, $p_lrt, $p_ct, $fdr_fcpt, $fdr_lrt, $fdr_ct) = split(/\t/);
	$smg_rel_pvalue{$gene} = $p_ct;
}
close(R);
INFO(scalar(keys(%smg_rel_pvalue))." genes read from file $smg_rel");

# read significantly mutated pathways at diagnosis
my %smp_dia_pvalue;
open(D,"$smp_dia") or croak "ERROR: could not read file $smp_dia\n";
<D>; # skip header
while(<D>)
{
	chomp;
	my ($id, $name, $class, $samples_affected, $total_variations, $p, $fdr, $num_genes, $genes) = split(/\t/);
	next if ($class ne "BBID" and $class ne "BIOCARTA" and $class ne "KEGG_PATHWAY" and $class ne "OMIM_DISEASE");
	next if ($p > 0.05);
	foreach my $g (split(",", $genes))
	{
		$smp_dia_pvalue{$g}{'pvalue'} = $p if (!exists $smp_dia_pvalue{$g}{'pvalue'} or $smp_dia_pvalue{$g}{'pvalue'} > $p);
		$smp_dia_pvalue{$g}{$id."|".$name} = $p;
	}
}
close(D);
INFO(scalar(keys(%smp_dia_pvalue))." pathways read from file $smp_dia");

# read significantly mutated pathways at relapse
my %smp_rel_pvalue;
open(D,"$smp_rel") or croak "ERROR: could not read file $smp_rel\n";
<D>; # skip header
while(<D>)
{
	chomp;
	my ($id, $name, $class, $samples_affected, $total_variations, $p, $fdr, $num_genes, $genes) = split(/\t/);
	next if ($class ne "BBID" and $class ne "BIOCARTA" and $class ne "KEGG_PATHWAY" and $class ne "OMIM_DISEASE");
	next if ($p > 0.05);
	foreach my $g (split(",", $genes))
	{
		$smp_rel_pvalue{$g}{'pvalue'} = $p if (!exists $smp_rel_pvalue{$g}{'pvalue'} or $smp_rel_pvalue{$g}{'pvalue'} > $p);
		$smp_rel_pvalue{$g}{$id."|".$name} = $p;
	}
}
close(D);
INFO(scalar(keys(%smp_rel_pvalue))." pathways read from file $smp_rel");

open(M,"$gene_patient_matrix") or croak "ERROR: could not read file $gene_patient_matrix\n";
my (%gene_info, @patients_dia, @patients_rel, @patients_cons);
my $header = <M>;
chomp($header);
my @hfields = split("\t", $header);
for (my $d = 9; $d <= 28; $d ++) { push(@patients_dia, $hfields[$d]); }
for (my $d = 33; $d <= 52; $d ++) { push(@patients_rel, $hfields[$d]); }
for (my $d = 55; $d <= 74; $d ++) { push(@patients_cons, $hfields[$d]); }

while(<M>)
{
	chomp;
	my @fields = split /\t/;

	my $g = $fields[0];

	$gene_info{$g}{'desc'} = $fields[1];
	$gene_info{$g}{'exons'} = $fields[2];
	$gene_info{$g}{'tr_len'} = $fields[3];
	$gene_info{$g}{'cds_len'} = $fields[4];
		
	$gene_info{$g}{'freq-dia'} = $fields[5];
	$gene_info{$g}{'tot-dia'} = $fields[6];
	$gene_info{$g}{'freq-dia-ns'} = $fields[7];
	$gene_info{$g}{'tot-dia-ns'} = $fields[8];
	for (my $d = 9; $d <= 28; $d ++) { $gene_info{$g}{$patients_dia[$d-9]} = $fields[$d]; }

	$gene_info{$g}{'freq-rel'} = $fields[29];
	$gene_info{$g}{'tot-rel'} = $fields[30];
	$gene_info{$g}{'freq-rel-ns'} = $fields[31];
	$gene_info{$g}{'tot-rel-ns'} = $fields[32];
	for (my $d = 33; $d <= 52; $d ++) { $gene_info{$g}{$patients_rel[$d-33]} = $fields[$d]; }

	$gene_info{$g}{'freq-cons'} = $fields[53];
	$gene_info{$g}{'tot-cons'} = $fields[54];
	for (my $d = 55; $d <= 74; $d ++) { $gene_info{$g}{$patients_cons[$d-55]} = $fields[$d]; }
}
close(M);
INFO(scalar(keys(%gene_info))." genes read from file $gene_patient_matrix");

# output header
print "gene\tdescr\texons\ttr_len\tcds_len\t";
print "freq-dia\ttot-dia\tfreq-dia-ns\ttot-dia-ns\tp-gene-dia\tp-pw-dia";
map { print "\t$_" } (@patients_dia);
print "\tfreq-rel\ttot-rel\tfreq-rel-ns\ttot-rel-ns\tp-gene-rel\tp-pw-rel";
map { print "\t$_" } (@patients_rel);
print "\tfreq-cons\ttot-cons";
map { print "\t$_" } (@patients_cons);
print "\tenr-pw-dia\tenr-pw-rel";
print "\n";

foreach my $g (keys(%gene_info))
{
	print "$g";
	print "\t",$gene_info{$g}{'desc'},"\t",$gene_info{$g}{'exons'},"\t",$gene_info{$g}{'tr_len'},"\t",$gene_info{$g}{'cds_len'};

	print "\t".(defined $gene_info{$g}{'freq-dia'} ? $gene_info{$g}{'freq-dia'} : "");
	print "\t".(defined $gene_info{$g}{'tot-dia'} ? $gene_info{$g}{'tot-dia'} : "");
	print "\t".(defined $gene_info{$g}{'freq-dia-ns'} ? $gene_info{$g}{'freq-dia-ns'} : "");
	print "\t".(defined $gene_info{$g}{'tot-dia-ns'} ? $gene_info{$g}{'tot-dia-ns'} : "");
	print "\t".(defined $smg_dia_pvalue{$g} ? $smg_dia_pvalue{$g} : ""); 
	print "\t".(defined $smp_dia_pvalue{$g}{'pvalue'} ? $smp_dia_pvalue{$g}{'pvalue'} : ""); 
	map { print "\t".(defined $gene_info{$g}{$_} ? $gene_info{$g}{$_} : "") } (@patients_dia);

	print "\t".(defined $gene_info{$g}{'freq-rel'} ? $gene_info{$g}{'freq-rel'} : "");
	print "\t".(defined $gene_info{$g}{'tot-rel'} ? $gene_info{$g}{'tot-rel'} : "");
	print "\t".(defined $gene_info{$g}{'freq-rel-ns'} ? $gene_info{$g}{'freq-rel-ns'} : "");
	print "\t".(defined $gene_info{$g}{'tot-rel-ns'} ? $gene_info{$g}{'tot-rel-ns'} : "");
	print "\t".(defined $smg_rel_pvalue{$g} ? $smg_rel_pvalue{$g} : ""); 
	print "\t".(defined $smp_rel_pvalue{$g}{'pvalue'} ? $smp_rel_pvalue{$g}{'pvalue'} : ""); 
	map { print "\t".(defined $gene_info{$g}{$_} ? $gene_info{$g}{$_} : "") } (@patients_rel);

	print "\t".(defined $gene_info{$g}{'freq-cons'} ? $gene_info{$g}{'freq-cons'} : "");
	print "\t".(defined $gene_info{$g}{'tot-cons'} ? $gene_info{$g}{'tot-cons'} : "");
	map { print "\t".(defined $gene_info{$g}{$_} ? $gene_info{$g}{$_} : "") } (@patients_cons);

	print "\t";
	my @enriched_pw_dia;
	map { push(@enriched_pw_dia, $_."(".$smp_dia_pvalue{$g}{$_}.")") if ($_ ne 'pvalue') } sort {$smp_dia_pvalue{$g}{$a} <=> $smp_dia_pvalue{$g}{$b}} keys(%{$smp_dia_pvalue{$g}});
	print join(",", @enriched_pw_dia);

	print "\t";
	my @enriched_pw_rel;
	map { push(@enriched_pw_rel, $_."(".$smp_rel_pvalue{$g}{$_}.")") if ($_ ne 'pvalue') } sort {$smp_rel_pvalue{$g}{$a} <=> $smp_rel_pvalue{$g}{$b}} keys(%{$smp_rel_pvalue{$g}});
	print join(",", @enriched_pw_rel);
	 
	print "\n";
}
