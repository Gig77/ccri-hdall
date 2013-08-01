use warnings FATAL => qw( all );
use strict;

use Carp;
use Getopt::Long;

# parse detailed results first
my ($sm_pathways, $sm_pathways_detail, $pathway_file, $gene_pathway_matrix);
GetOptions
(
	"pathway-file=s" => \$pathway_file, # input pathway file for music path-scan
	"sm-pathways=s" => \$sm_pathways, # result list from genome music path-scan  
	"sm-pathways-detail=s" => \$sm_pathways_detail, # detaild result from genome music path-scan (= second output file)  
	"gene-pathway-matrix=s" => \$gene_pathway_matrix # OUTPUT: gene-pathway-matrix  
);

# read pathway input file
my %pathways;
open(P, "$pathway_file") or die "ERROR: could not open pathway file $pathway_file\n";
while(<P>)
{
	chomp;
	my ($path_id, $path_name, $class, $gene_line, $diseases, $drugs, $description) = split(/\t/);
	$pathways{$path_id}{'size'} = scalar(split('\|', $gene_line));
}
close(P);

my $content = "";
open(D,"$sm_pathways_detail") or die "ERROR: could not read file $sm_pathways_detail\n";
while(<D>) 
{
#	chomp;
	$content .= $_;
}
close(D);

my (%pw, %genes);
while ($content =~ /Pathway: (.*?)\nName: (.*?)\nClass: (.*?)\nP-value: (.*?)\nFDR: (.*?)\nDescription:(.*?)\n(.*?)\nSamples with mutations \(#hits\): (.*?)\n/sg)
{
#	next if ($1 ne 'GOTERM_BP_FAT:GO:0008284');
#	print "Pathway: $1\nName: $2\nClass: $3\nP-value: $4\nFDR: $5\nDescription: $6\n";
		
	my $id = $1;
	$pw{$id}{name} = $2;
	$pw{$id}{class} = $3;
	$pw{$id}{pvalue} = $4;
	$pw{$id}{fdr} = $5;
	$pw{$id}{desc} = $6;
	
	my $genes_per_sample = $7;
	foreach my $g (split(/\n/, $genes_per_sample))
	{
		my ($sample, $genes) = split (":", $g);
#		print "  $g\n";
		$pw{$id}{genes_per_sample}{$sample} = $genes;
		
		foreach my $gene (split(",", $genes))
		{
			$pw{$id}{genes}{$gene} = $pw{$id}{genes}{$gene} ? $pw{$id}{genes}{$gene} + 1 : 1;
			$genes{$gene}{$id} = 1;		
		}
	}
}	

# TABLE: sm_pathways.annotated
print "Pathway\tName\tClass\tSize\tSamples_Affected\tTotal_Variations\tp-value\tFDR\tNumGenes\tGenes\n";
open(D,"$sm_pathways") or die "ERROR: could not read file $sm_pathways\n";
<D>; # skip header
while(<D>) 
{
	chomp;
	my ($pathway, $name, $class, $num_samples, $num_genes, $pvalue, $fdr) = split(/\t/);
	print "$pathway\t$name\t$class\t",$pathways{$pathway}{'size'},"\t$num_samples\t$num_genes\t$pvalue\t$fdr\t";
	print scalar(keys(%{$pw{$pathway}{genes}})),"\t";
	my @genes_sorted = sort { $pw{$pathway}{genes}{$b} <=> $pw{$pathway}{genes}{$a} } keys(%{$pw{$pathway}{genes}});
	my @genes_print = map { "$_(".$pw{$pathway}{genes}{$_}.")"  } @genes_sorted; 
	print join(",", @genes_print),"\n";
}
close(D);

if ($gene_pathway_matrix)
{
	open(M,">$gene_pathway_matrix") or croak("ERROR: could not write to file $gene_pathway_matrix");
	map { print M "\t$_" } keys(%pw);
	print M "\n";
	foreach my $g (keys(%genes))
	{
		print M "$g";
		map { print M $genes{$g}{$_} ? "\t1" : "\t0" } keys(%pw);
		print M "\n";
	}	
	close(M);
}
