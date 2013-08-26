use warnings FATAL => qw( all );
use strict;

use Carp;
use Getopt::Long;
use List::Util qw(min max);

# STDIN: list with enriched pathways from DAVID
# STDOUT: pathway/patient matrix indicating the number of times a pathway is mutated across patients

my ($max_pvalue, $max_pathway_size, $subtract_pathway_file);
GetOptions
(
	"max-p-value=f" => \$max_pvalue,  
	"max-pathway-size=i" => \$max_pathway_size,
	"subtract-pathways=s" => \$subtract_pathway_file  
);

croak "ERROR: --max-p-value not specified" if (!defined $max_pvalue);
croak "ERROR: --max-pathway-size not specified" if (!defined $max_pathway_size);

my %subtract_pathways;
my $subtract_pvalue = 0.01;
if ($subtract_pathway_file)
{
	# TABLE: sm_pathways.annotated
	open(S, "$subtract_pathway_file") or croak "ERROR: could not read file $subtract_pathway_file";
	<S>; # skip header
	while(<S>)
	{
		my ($pathway, $name, $class, $size, $samples, $tot_var, $p_value, $fdr, $num_genes, $genes) = split /\t/;
		$subtract_pathways{$pathway} = 1 if ($p_value < $subtract_pvalue);
	}
	close(S);
}
print STDERR "INFO: Will exclude ".keys(%subtract_pathways)." pathways with p-value < $subtract_pvalue in diagnosis\n";

# TABLE: sm_pathways.annotated
my (%pw, %genes);
<>; # skip header
while(<>)
{
	chomp;
	my ($pathway, $name, $class, $size, $samples, $tot_var, $p_value, $fdr, $num_genes, $genes) = split /\t/;
	croak "ERROR: Could not parse pathway entry: $_\n" if (!$pathway or !$name or !$class or !$genes);

	next if ($p_value > $max_pvalue);
	next if ($size > $max_pathway_size);
#	next unless ($class =~ /^c\d\_/); # only MSigDB pathways
	next if ($class =~ /^c1_positional/);
#	next if ($name =~ /^(REACTOME|KEGG|BIOCARTA|PID)_/);
	next if ($subtract_pathways{$pathway});
	
	$pw{$pathway}{name} = $name;
	$pw{$pathway}{class} = $class;
	$pw{$pathway}{size} = $size;
	$pw{$pathway}{pvalue} = $p_value;

	foreach my $g (split(",", $genes))
	{
		$g =~ s/\(\d+\)//;
		$genes{$g}{$pathway} = 1;
	}
}

map { print "\t",$pw{$_}{class},"|",$pw{$_}{name},"|",$pw{$_}{size},"|",sprintf("%.2g", $pw{$_}{pvalue}) } keys(%pw);
print "\n";
foreach my $g (keys(%genes))
{
	print "$g";
	map { print $genes{$g}{$_} ? "\t1" : "\t0" } keys(%pw);
	print "\n";
}	
