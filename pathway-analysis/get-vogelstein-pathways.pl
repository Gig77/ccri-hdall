use warnings FATAL => qw( all );
use strict;

my %pathways;

# read pathways
open(G,"/mnt/projects/hdall/data/vogelstein_2013/1235122TablesS1-4.tsv") or die "could not open file /mnt/projects/hdall/data/vogelstein_2013/1235122TablesS1-4.tsv";
<G>; # skip header
while(<G>)
{
	chomp;
	my ($symbol, $name, $num_mut_tum, $oncogene_score, $suppressor_score, $classification, $core_pathway, $process) = split(/\t/);
	$pathways{$symbol} = $core_pathway;
}
close(G);

while(<>)
{
	chomp;
	my ($gene, $status) = split("\t");
	next if ($gene eq "gene");
	
	print "$gene\t$status\t";
	print exists $pathways{$gene} ? $pathways{$gene} : "";
	print "\n"; 
}