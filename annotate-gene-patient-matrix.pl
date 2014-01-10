use warnings FATAL => qw( all );
use strict;

use lib "$ENV{HOME}/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Carp;
use Getopt::Long;

# parse detailed results first
my ($gene_patient_matrix, $smg_dia, $smg_rel, $smp_dia_file, $smp_rel_file, $cnv, $gene_pathway_matrix, $gene_panel_file, $max_pvalue, $max_pathway_size, $curated_freq_file);
GetOptions
(
	"gene-patient-matrix=s" => \$gene_patient_matrix,
	"smg-dia=s" => \$smg_dia,   
	"smg-rel=s" => \$smg_rel,   
	"smp-dia=s" => \$smp_dia_file,  
	"smp-rel=s" => \$smp_rel_file,
	"cnv=s" => \$cnv,
	"gene-pathway-matrix=s" => \$gene_pathway_matrix, # clustered gene-pathway-matrix to determine order of genes and pathways
	"panel-genes=s" => \$gene_panel_file, # names of genes selected for resequencing   
	"max-p-value=f" => \$max_pvalue,  
	"max-pathway-size=i" => \$max_pathway_size,
	"curated-frequencies=s" => \$curated_freq_file
);

croak "ERROR: Pathway group file not specified" if (!$gene_pathway_matrix);

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
croak "ERROR: no significanly mutated genes found in file $smg_rel" if (keys(%smg_rel_pvalue) == 0);

# read significantly mutated pathways at diagnosis
# TABLE: sm_pathways.annotated
my (%smp_dia, %smp_dia_genes);
open(D,"$smp_dia_file") or croak "ERROR: could not read file $smp_dia_file\n";
<D>; # skip header
while(<D>)
{
	chomp;
	my ($id, $name, $class, $size, $samples_affected, $total_variations, $p, $fdr, $num_genes, $genes, $link) = split(/\t/);
	next if ($size > $max_pathway_size);
#	next if ($class ne "NCI");
	next if ($p > $max_pvalue);
	next if ($class =~ /(c1_|c3_|c7_|c4_)/);

	$smp_dia{"$class|$name|$size"} = $p;
	foreach my $g (split(",", $genes))
	{
		$g =~ s/\(\d+\)//;
		if (!exists $smp_dia_genes{$g}{'pvalue'} or $smp_dia_genes{$g}{'pvalue'} > $p)
		{
			$smp_dia_genes{$g}{'pvalue'} = $p;
			$smp_dia_genes{$g}{'pathway'} = "$class|$name|$size";
		}
		$smp_dia_genes{$g}{"$class|$name|$size"} = $p;
	}
}
close(D);
INFO(scalar(keys(%smp_dia_genes))." pathways read from file $smp_dia_file");
INFO("WARNING: no significanly mutated pathways found in file $smp_dia_file") if (keys(%smp_dia_genes) == 0);

# read significantly mutated pathways at relapse
# TABLE: sm_pathways.annotated
my (%smp_rel, %smp_rel_genes);
open(D,"$smp_rel_file") or croak "ERROR: could not read file $smp_rel_file\n";
<D>; # skip header
while(<D>)
{
	chomp;
	my ($id, $name, $class, $size, $samples_affected, $total_variations, $p, $fdr, $num_genes, $genes, $link) = split(/\t/);
	next if ($size > $max_pathway_size);
#	next if ($class ne "NCI");
	next if ($p > $max_pvalue);
	next if ($class =~ /(c1_|c3_|c7_|c4_)/);

	$smp_rel{"$class|$name|$size"} = $p;
	foreach my $g (split(",", $genes))
	{
		$g =~ s/\(\d+\)//;
		if (!exists $smp_rel_genes{$g}{'pvalue'} or $smp_rel_genes{$g}{'pvalue'} > $p)
		{
			$smp_rel_genes{$g}{'pvalue'} = $p;
			$smp_rel_genes{$g}{'pathway'} = "$class|$name|$size";
		}
		$smp_rel_genes{$g}{"$class|$name|$size"} = $p;
	}
}
close(D);
INFO(scalar(keys(%smp_rel_genes))." pathways read from file $smp_rel_file");
INFO("WARNING: no significanly mutated pathways found in file $smp_rel_file") if (keys(%smp_rel_genes) == 0);

# read copy-number info
my (%cnvs);
my $cnv_read = 0;
open(D,"$cnv") or croak "ERROR: could not read CNV file $cnv\n";
<D>; # skip header: patient\tsample\tgene\tchr\tstart\tend\ttrlen\tcdslen\texons\tcosmic\tdescr\tcnumber\tevent\tevent_coordinate\tevent_size\tnum_genes\n
while(<D>)
{
	chomp;
	my ($patient, $sample, $gene, $chr, $start, $end, $trlen, $cdslen, $exons, $cosmic, $descr, $cnumber, $event, $event_coordinate, $event_size, $num_genes) = split(/\t/);
	$cnvs{$patient}{$sample}{$gene}{'cnumber'} = $cnumber;
	$cnvs{$patient}{$sample}{$gene}{'event'} = $event;
	$cnv_read ++;
}
close(D);
INFO("$cnv_read CNVs read from file $cnv");

# read gene/pathway matrix
my %gpm;
my %funct_gene_order;
my $lines_read = 0;
open(PG, "$gene_pathway_matrix") or croak "ERROR: could not open file $gene_pathway_matrix\n"; 
my $hline = <PG>; 
chomp($hline);
my @pw_names = split("\t", $hline);
shift(@pw_names);
for (my $i = 0; $i < @pw_names; $i ++) { my ($pwcategory, $pwname, $pwsize) = split ('\|', $pw_names[$i]); $pw_names[$i] = "$pwcategory|$pwname|$pwsize"; }
while(<PG>)
{
	chomp;
	my ($gene, @cells) = split /\t/;
	croak "ERROR: coud not parse following line:\n$_\n" if (@cells < @pw_names);
	
	$funct_gene_order{$gene} = $lines_read+1;
	for (my $i = 0; $i < @cells; $i ++)
	{
		$gpm{$gene}{$pw_names[$i]} = 1 if ($cells[$i]);
	}
	$lines_read ++;
}
close(PG);
INFO("$lines_read lines read from file $gene_pathway_matrix");

# read panel genes
my %panel_genes;
if ($gene_panel_file)
{
	open(P,"$gene_panel_file") or croak "ERROR: could not read file $gene_panel_file\n";
	while(<P>)
	{
		chomp;
		$panel_genes{$_} = 1;
	}
	close(P);
	INFO(scalar(keys(%panel_genes))." genes read from file $gene_panel_file");	
}

# read curated frequency number
# TABLE: freq-rel-ns-curated.tsv
my %curated_freq;
if ($curated_freq_file)
{
	open(C,"$curated_freq_file") or croak "ERROR: could not read CNV file $curated_freq_file\n";
	<C>; # skip header:
	while(<C>)
	{
		chomp;
		my ($gene, $freq, $comment) = split(/\t/);
		$curated_freq{$gene} = $freq;
	}
	close(C);
	INFO("".scalar(keys(%curated_freq))." genes read from file $curated_freq_file");	
}

# TABLE: gene-patient-matrix
open(M,"$gene_patient_matrix") or croak "ERROR: could not read file $gene_patient_matrix\n";
my (%gene_info, @patients_dia, @patients_rel, @patients_cons);
my $header = <M>;
chomp($header);
my @hfields = split("\t", $header);
for (my $d = 18; $d <= 37; $d ++) { push(@patients_dia, $hfields[$d]); }
for (my $d = 47; $d <= 66; $d ++) { push(@patients_rel, $hfields[$d]); }
for (my $d = 70; $d <= 89; $d ++) { push(@patients_cons, $hfields[$d]); }

while(<M>)
{
	chomp;
	my @fields = split /\t/;

	my $g = $fields[0];

	$gene_info{$g}{'desc'} = $fields[1];
	$gene_info{$g}{'chr'} = $fields[2];
	$gene_info{$g}{'start'} = $fields[3];
	$gene_info{$g}{'end'} = $fields[4];
	$gene_info{$g}{'exons'} = $fields[5];
	$gene_info{$g}{'tr_len'} = $fields[6];
	$gene_info{$g}{'cds_len'} = $fields[7];
	$gene_info{$g}{'cosmic'} = $fields[8];
		
	$gene_info{$g}{'freq-dia'} = $fields[9];
	$gene_info{$g}{'tot-dia'} = $fields[10];
	$gene_info{$g}{'freq-dia-ns'} = $fields[11];
	$gene_info{$g}{'freq-dia-ns-af20'} = $fields[12];
	$gene_info{$g}{'tot-dia-ns'} = $fields[13];
	$gene_info{$g}{'max-af-dia'} = $fields[14];
	$gene_info{$g}{'max-af-dia-ns'} = $fields[15];
	$gene_info{$g}{'imp-ex-dia'} = $fields[16];
	$gene_info{$g}{'imp-ex-dia-ns'} = $fields[17];
	for (my $d = 18; $d <= 37; $d ++) { $gene_info{$g}{$patients_dia[$d-18]} = $fields[$d]; }

	$gene_info{$g}{'freq-rel'} = $fields[38];
	$gene_info{$g}{'tot-rel'} = $fields[39];
	$gene_info{$g}{'freq-rel-ns'} = $fields[40];
	$gene_info{$g}{'freq-rel-ns-af20'} = $fields[41];
	$gene_info{$g}{'tot-rel-ns'} = $fields[42];
	$gene_info{$g}{'max-af-rel'} = $fields[43];
	$gene_info{$g}{'max-af-rel-ns'} = $fields[44];
	$gene_info{$g}{'imp-ex-rel'} = $fields[45];
	$gene_info{$g}{'imp-ex-rel-ns'} = $fields[46];
	for (my $d = 47; $d <= 66; $d ++) { $gene_info{$g}{$patients_rel[$d-47]} = $fields[$d]; }

	$gene_info{$g}{'freq-cons'} = $fields[67];
	$gene_info{$g}{'freq-cons-raise'} = $fields[68];
	$gene_info{$g}{'tot-cons'} = $fields[69];
	for (my $d = 70; $d <= 89; $d ++) { $gene_info{$g}{$patients_cons[$d-70]} = $fields[$d]; }
	
	$gene_info{$g}{'imp-domains-dia'} = $fields[90];
	$gene_info{$g}{'imp-domains-rel'} = $fields[91];
}
close(M);
INFO(scalar(keys(%gene_info))." genes read from file $gene_patient_matrix");

# TABLE: gene-patient-matrix.annotated

# output header
print "gene\ton-panel\tdescr\tchr\tstart\tend\texons\ttr_len\tcds_len\tcosmic\t";
print "freq-dia\ttot-dia\tfreq-dia-ns\tfreq-dia-ns-af20\ttot-dia-ns\tmax-af-dia\tmax-af-dia-ns\timp-ex-dia\timp-ex-dia-ns\tp-gene-dia\tp-pw-dia\ttop-pw-dia";
map { print "\t$_" } (@patients_dia);
print "\tfreq-rel\ttot-rel\tfreq-rel-ns\tfreq-rel-ns-af20\tfreq-rel-ns-curated\ttot-rel-ns\tmax-af-rel\tmax-af-rel-ns\timp-ex-rel\timp-ex-rel-ns\tp-gene-rel\tp-pw-rel\ttop-pw-rel";
map { print "\t$_" } (@patients_rel);
print "\tfreq-cons\tfreq-cons-raise\ttot-cons";
map { print "\t$_" } (@patients_cons);
print "\timp-domains-dia\timp-domains-rel";
print "\torder";

my @pw_names_rel;
map { if ($smp_rel{$_}) { print "\t$_|",sprintf("%.1e", $smp_rel{$_}); push(@pw_names_rel, $_); } } (@pw_names);
print "\n";

foreach my $g (keys(%gene_info))
{
	print "$g";
	print "\t", $panel_genes{$g} ? "yes" : "-";
	print "\t",$gene_info{$g}{'desc'},"\t",$gene_info{$g}{'chr'},"\t",$gene_info{$g}{'start'},"\t",$gene_info{$g}{'end'},"\t",$gene_info{$g}{'exons'},"\t",$gene_info{$g}{'tr_len'},"\t",$gene_info{$g}{'cds_len'},"\t",$gene_info{$g}{'cosmic'};

	print "\t".(defined $gene_info{$g}{'freq-dia'} ? $gene_info{$g}{'freq-dia'} : "");
	print "\t".(defined $gene_info{$g}{'tot-dia'} ? $gene_info{$g}{'tot-dia'} : "");
	print "\t".(defined $gene_info{$g}{'freq-dia-ns'} ? $gene_info{$g}{'freq-dia-ns'} : "");
	print "\t".(defined $gene_info{$g}{'freq-dia-ns-af20'} ? $gene_info{$g}{'freq-dia-ns-af20'} : "");
	print "\t".(defined $gene_info{$g}{'tot-dia-ns'} ? $gene_info{$g}{'tot-dia-ns'} : "");
	print "\t".(defined $gene_info{$g}{'max-af-dia'} ? $gene_info{$g}{'max-af-dia'} : "");
	print "\t".(defined $gene_info{$g}{'max-af-dia-ns'} ? $gene_info{$g}{'max-af-dia-ns'} : "");
	print "\t".(defined $gene_info{$g}{'imp-ex-dia'} ? $gene_info{$g}{'imp-ex-dia'} : "");
	print "\t".(defined $gene_info{$g}{'imp-ex-dia-ns'} ? $gene_info{$g}{'imp-ex-dia-ns'} : "");

	print "\t".(defined $smg_dia_pvalue{$g} ? $smg_dia_pvalue{$g} : ""); 
	print "\t".(defined $smp_dia_genes{$g}{'pvalue'} ? $smp_dia_genes{$g}{'pvalue'} : "");
	print "\t".(defined $smp_dia_genes{$g}{'pathway'} ? $smp_dia_genes{$g}{'pathway'} : ""); 

	foreach my $pdia (@patients_dia)
	{
		print "\t";
		next if (!$gene_info{$g}{$pdia} or $gene_info{$g}{$pdia} eq " ");
		print $gene_info{$g}{$pdia};
		
		my ($p) = $pdia =~ /(.*?)-dia/; 
		if ($cnvs{$p}{'rem_dia'}{$g})
		{
			print "|".$cnvs{$p}{'rem_dia'}{$g}{'event'}.":".$cnvs{$p}{'rem_dia'}{$g}{'cnumber'};
		}		
	} 

	print "\t".(defined $gene_info{$g}{'freq-rel'} ? $gene_info{$g}{'freq-rel'} : "");
	print "\t".(defined $gene_info{$g}{'tot-rel'} ? $gene_info{$g}{'tot-rel'} : "");
	print "\t".(defined $gene_info{$g}{'freq-rel-ns'} ? $gene_info{$g}{'freq-rel-ns'} : "");
	print "\t".(defined $gene_info{$g}{'freq-rel-ns-20'} ? $gene_info{$g}{'freq-rel-ns-20'} : "");
	print "\t".(defined $curated_freq{$g} ? $curated_freq{$g} : (defined $gene_info{$g}{'freq-rel-ns'} ? $gene_info{$g}{'freq-rel-ns'} : ""));
	print "\t".(defined $gene_info{$g}{'tot-rel-ns'} ? $gene_info{$g}{'tot-rel-ns'} : "");
	print "\t".(defined $gene_info{$g}{'max-af-rel'} ? $gene_info{$g}{'max-af-rel'} : "");
	print "\t".(defined $gene_info{$g}{'max-af-rel-ns'} ? $gene_info{$g}{'max-af-rel-ns'} : "");
	print "\t".(defined $gene_info{$g}{'imp-ex-rel'} ? $gene_info{$g}{'imp-ex-rel'} : "");
	print "\t".(defined $gene_info{$g}{'imp-ex-rel-ns'} ? $gene_info{$g}{'imp-ex-rel-ns'} : "");

	print "\t".(defined $smg_rel_pvalue{$g} ? $smg_rel_pvalue{$g} : ""); 
	print "\t".(defined $smp_rel_genes{$g}{'pvalue'} ? $smp_rel_genes{$g}{'pvalue'} : ""); 
	print "\t".(defined $smp_rel_genes{$g}{'pathway'} ? $smp_rel_genes{$g}{'pathway'} : ""); 

	foreach my $prel (@patients_rel)
	{
		print "\t";
		next if (!$gene_info{$g}{$prel} or $gene_info{$g}{$prel} eq " ");
		print $gene_info{$g}{$prel};
		
		my ($p) = $prel =~ /(.*?)-rel/; 
		if ($cnvs{$p}{'rem_rel'}{$g})
		{
			print "|".$cnvs{$p}{'rem_rel'}{$g}{'event'}.":".$cnvs{$p}{'rem_rel'}{$g}{'cnumber'};
		}		
	} 

	print "\t".(defined $gene_info{$g}{'freq-cons'} ? $gene_info{$g}{'freq-cons'} : "");
	print "\t".(defined $gene_info{$g}{'freq-cons-raise'} ? $gene_info{$g}{'freq-cons-raise'} : "");
	print "\t".(defined $gene_info{$g}{'tot-cons'} ? $gene_info{$g}{'tot-cons'} : "");
	map { print "\t".(defined $gene_info{$g}{$_} ? $gene_info{$g}{$_} : "") } (@patients_cons);

	print "\t".(defined $gene_info{$g}{'imp-domains-dia'} ? $gene_info{$g}{'imp-domains-dia'} : "");
	print "\t".(defined $gene_info{$g}{'imp-domains-rel'} ? $gene_info{$g}{'imp-domains-rel'} : "");

	print "\t".(defined $funct_gene_order{$g} ? $funct_gene_order{$g} : "");

	map { print "\t", defined $smp_rel_genes{$g}{$_} ? "x" : "" } (@pw_names_rel);
		
#	print "\t";
#	my @enriched_pw_dia;
#	map { push(@enriched_pw_dia, $_."(".sprintf("%.1e", $smp_dia_genes{$g}{$_}).")") if ($_ ne 'pvalue') } sort {$smp_dia_genes{$g}{$a} <=> $smp_dia_genes{$g}{$b}} keys(%{$smp_dia_genes{$g}});
#	print join(",", @enriched_pw_dia);
#
#	print "\t";
#	my @enriched_pw_rel;
#	map { push(@enriched_pw_rel, $_."(".sprintf("%.1e", $smp_rel_genes{$g}{$_}).")") if ($_ ne 'pvalue') } sort {$smp_rel_genes{$g}{$a} <=> $smp_rel_genes{$g}{$b}} keys(%{$smp_rel_genes{$g}});
#	print join(",", @enriched_pw_rel);

#	print "\t";
#	my @enriched_pw_rel_spec;
#	map { push(@enriched_pw_rel_spec, $_."(".sprintf("%.1e", $smp_rel{$_}).")") if ($_ ne 'pvalue' and $smp_rel{$_} < 0.01 and (!defined $smp_dia{$_} or $smp_dia{$_} > 0.1)) } sort {$smp_rel_genes{$g}{$a} <=> $smp_rel_genes{$g}{$b}} keys(%{$smp_rel_genes{$g}});
#	print join(",", @enriched_pw_rel_spec);
	 
	print "\n";
}
