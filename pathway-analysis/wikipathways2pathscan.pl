use warnings FATAL => qw( all );
use strict;
use Carp;

use lib "$ENV{HOME}/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Getopt::Long;

my ($wikidir);
GetOptions
(
	"wiki-dir=s" => \$wikidir  # directory containing WikiPathways (one txt file for each pathway)
);
croak "ERROR: --wiki-dir not specified" if (!$wikidir);

# read id mapping
my %id2sym;
open(M, "$ENV{HOME}/hdall/results/id-mappings.tsv") or croak "ERROR: could not read id mappings\n";
while(<M>)
{
	chomp;
	my ($sym, $id) = split(/\t/);
	$id2sym{$id} = $sym;
}
close(M);
INFO(scalar(keys(%id2sym))." id mappgins read from file $ENV{HOME}/hdall/results/id-mappings.tsv");

# read biomart id mapping to get entrez ids and additional uniprot ids
my %sym2entrez;
open(G, "$ENV{HOME}/hdall/results/gene-id-mapping.biomart-0.7.tsv") or die "ERROR: could not read gene list\n";
while(<G>)
{
	chomp;
	my ($approved_symbol, $entrez_gene_id, $accession_numbers, $approved_name, 
		$previous_symbols, $previous_names, $aliases, $name_aliases, $entrez_gene_id_ncbi, $ensembl_id, $uniprot, $refseq, $ucsc) = split("\t");

	$id2sym{$entrez_gene_id_ncbi} = $approved_symbol if ($entrez_gene_id_ncbi);
	$id2sym{$entrez_gene_id} = $approved_symbol if ($entrez_gene_id);
	$id2sym{$ensembl_id} = $approved_symbol if ($ensembl_id);
}
close(G);

#	# MuSiC required format:
#	# This is a tab-delimited file prepared from a pathway database (such as KEGG), with the columns: 
#	# [path_id, path_name, class, gene_line, diseases, drugs, description] The latter three columns 
#	# are optional (but are available on KEGG). The gene_line contains the "entrez_id:gene_name" of 
#	# all genes involved in this pathway, each separated by a "|" symbol.
#	# Example:
#	# hsa00061      Fatty acid biosynthesis	      Lipid Metabolism     31:ACACA|32:ACACB|27349:MCAT|2194:FASN|54995:OXSM|55301:OLAH	

opendir(D, $wikidir) or croak "ERROR: could not read directory $wikidir!";
while (my $f = readdir(D))
{
	if ($f =~ /Hs_(.*?)\.txt$/)
	{
		my $pw = $1;
		print "$pw\t$pw\tWikiPathways\t";
		
		my %genes;
		open(P, "$wikidir/$f") or croak "ERROR: could not read file $f";
		<P>;
		while (<P>)
		{
			chomp;
			my ($id, $db) = split("\t");
			
			my $sym = $id2sym{$id};
			$sym = $id2sym{$sym} if ($sym and $id2sym{$sym});
			$genes{"0000:$sym"} = 1 if ($sym);
		}
		close(P);
		
		print join('|', keys(%genes));
		print "\t\t\t\n";
	}
}
close(D);
#my ($read, $not_mapped) = (0, 0);
#my %pathways;
#while(<>)
#{
#	chomp;
#	my ($uniprot, $name_long, $name_short, $id) = split(/\t/);
#	
#	my $symbol = $uni2sym{$uniprot};
#	if (!$symbol)
#	{
#		$not_mapped ++;
#		#print "ERROR: Could not map uniprot accession $uniprot:  \n$_\n";
#		next;
#	}
#	
#	$pathways{$name_short}{'long'} = $name_long;
#	$pathways{$name_short}{'id'} = $id;
#	$pathways{$name_short}{'genes'}{$uniprot} = $symbol;
#	$read ++;
#}
#WARN("Could not map $not_mapped out of $read UniProt IDs.\n") if ($not_mapped > 0);
#
#foreach my $p (keys(%pathways))
#{
#	# MuSiC required format:
#	# This is a tab-delimited file prepared from a pathway database (such as KEGG), with the columns: 
#	# [path_id, path_name, class, gene_line, diseases, drugs, description] The latter three columns 
#	# are optional (but are available on KEGG). The gene_line contains the "entrez_id:gene_name" of 
#	# all genes involved in this pathway, each separated by a "|" symbol.
#	# Example:
#	# hsa00061      Fatty acid biosynthesis	      Lipid Metabolism     31:ACACA|32:ACACB|27349:MCAT|2194:FASN|54995:OXSM|55301:OLAH	
#	print $pathways{$p}{'id'}, "\t";
#	print "$p:", $pathways{$p}{'long'}, "\t";
#	print "NCI\t";
#
#	my @genes = values(%{$pathways{$p}{'genes'}});
#	for (my $i = 0; $i < @genes; $i++)
#	{
#		my $entrezid = $sym2entrez{$genes[$i]};
#		$entrezid = "0000" if (!$entrezid);
#		
#		print "$entrezid:$genes[$i]";
#		print "|" if ($i < @genes-1);
#	}
#		
#	print "\t\t\t\n";
#}