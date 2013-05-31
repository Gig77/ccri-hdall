use warnings;
use strict;

use Carp;
use Getopt::Long;

# parse detailed results first
my ($sm_pathways, $sm_pathways_detail);
GetOptions
(
	"sm-pathways=s" => \$sm_pathways, # result list from genome music path-scan  
	"sm-pathways-detail=s" => \$sm_pathways_detail # detaild result from genome music path-scan (= second output file)  
);

my $content = "";
open(D,"$sm_pathways_detail") or die "ERROR: could not read file $sm_pathways_detail\n";
while(<D>) 
{
#	chomp;
	$content .= $_;
}
close(D);

my %data;
while ($content =~ /Pathway: (.*?)\nName: (.*?)\nClass: (.*?)\nP-value: (.*?)\nFDR: (.*?)\nDescription:(.*?)\n(.*?)\nSamples with mutations \(#hits\): (.*?)\n/sg)
{
#	next if ($1 ne 'GOTERM_BP_FAT:GO:0008284');
#	print "Pathway: $1\nName: $2\nClass: $3\nP-value: $4\nFDR: $5\nDescription: $6\n";
		
	my $id = $1;
	$data{$id}{name} = $2;
	$data{$id}{class} = $3;
	$data{$id}{pvalue} = $4;
	$data{$id}{fdr} = $5;
	$data{$id}{desc} = $6;
	
	my $genes_per_sample = $7;
	foreach my $g (split(/\n/, $genes_per_sample))
	{
		my ($sample, $genes) = split (":", $g);
#		print "  $g\n";
		$data{$id}{genes_per_sample}{$sample} = $genes;
		
		foreach my $gene (split(",", $genes))
		{
			$data{$id}{genes}{$gene} = 1;		
		}
	}
}	

open(D,"$sm_pathways") or die "ERROR: could not read file $sm_pathways\n";
my $header = <D>; chomp($header);
$header .= "\tNumGenes\tGenes";
print "$header\n";
while(<D>) 
{
	chomp;
	my ($pathway, $name, $class, $num_samples, $num_genes, $pvalue, $fdr) = split(/\t/);
	print "$pathway\t$name\t$class\t$num_samples\t$num_genes\t$pvalue\t$fdr\t";
	print scalar(keys(%{$data{$pathway}{genes}})),"\t";
	print join(",", sort(keys(%{$data{$pathway}{genes}}))),"\n";
}
close(D);
