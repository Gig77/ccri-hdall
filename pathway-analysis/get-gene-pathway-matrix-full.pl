use warnings FATAL => qw( all );
use strict;

use Carp;
use Getopt::Long;
use List::Util qw(min max);

# reads pathscan pathway file and writes gene-pathway-matrix for clustering

my ($pathway_file);
GetOptions
(
	"pathway-file=s" => \$pathway_file  
);

#	# MuSiC required format:
#	# This is a tab-delimited file prepared from a pathway database (such as KEGG), with the columns: 
#	# [path_id, path_name, class, gene_line, diseases, drugs, description] The latter three columns 
#	# are optional (but are available on KEGG). The gene_line contains the "entrez_id:gene_name" of 
#	# all genes involved in this pathway, each separated by a "|" symbol.
#	# Example:
#	# hsa00061      Fatty acid biosynthesis	      Lipid Metabolism     31:ACACA|32:ACACB|27349:MCAT|2194:FASN|54995:OXSM|55301:OLAH
my (%p2g, %g2p);	
open(P, "$pathway_file") or croak("ERROR: Could not open pathway file $pathway_file");
while(<P>)
{
	chomp;
	my ($id, $pathway, $category, $genes_str) = split("\t");
	foreach my $entry (split('\|', $genes_str))
	{
		my ($entrezid, $gene) = split(":", $entry);
		$p2g{"$id|$pathway|$category"}{$gene} = 1;
		$g2p{$gene}{"$id|$pathway|$category"} = 1;
	}
}
print STDERR "".scalar(keys(%g2p))." genes in ".scalar(keys(%p2g))." pathways read from file $pathway_file\n";
close(P);

map { print "\t$_" } keys(%p2g);
print "\n";
foreach my $g (keys(%g2p))
{
	print "$g";
	map { print $p2g{$_}{$g} ? "\t1" : "\t0" } keys(%p2g);
	print "\n";
}