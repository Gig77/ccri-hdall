use warnings FATAL => qw( all );
use strict;
use Carp;

use lib "/mnt/projects/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Getopt::Long;
use Data::Dumper;

# INPUT FORMAT
# NUCLEAR_CHROMOSOME	http://www.broadinstitute.org/gsea/msigdb/cards/NUCLEAR_CHROMOSOME	MKI67IP	CBX1	TERF2IP	TTN	CBX5	MCM7	STAG3	H2AFY	TUBG1	NOL6	UBE2I	MCM3	PURB	PURA	KLHDC3	RAD51	CHMP1A	
# SMARCE1	ZNF238	TIMELESS	ZMIZ2	JUN	FOXC1	DMC1	LRPPRC	PAM	HMGB2	H1FNT	TIPIN	POLA1	NUFIP1	CHEK1	RCC1	RPA4	RPA3	RPA1	RPA2	ACD	RGS12	NPM2	ERCC4	ERCC1	PIF1	SUV39H1	SYCE2	
# SYCE1	SMC2	LEPREL4	ATRX	ZBED1	H2AFY2	SMC1A	HDAC8	REPIN1

my ($entrez_id_file, $category);
GetOptions
(
	"entrez-id-file=s" => \$entrez_id_file,  # Current MSigDB gene sets, Entrez IDs
	"category=s" => \$category # MSigDB category
);

die "ERROR: --entrez-id-file not specified\n" if (!$entrez_id_file);
die "ERROR: Invalid --category\n" if (!$category or $category !~ /(c1_positional|c2_curated|c3_motif|c4_computational|c5_GO|c6_oncogenic_signature|c7_immunologic_signature)/);

# read and remember entrez-ids
open(E, "$entrez_id_file") or die "ERROR: Could not read entrez id file $entrez_id_file\n";
my %entrez_ids;
while(<E>)
{
	chomp;
	my @fields = split /\t/;
	my $pathway = shift(@fields);
	$pathway =~ s/\)//g;	
	my $url = shift(@fields);
	my @ids = @fields;
	
	for (my $i = 0; $i < @ids; $i++)
	{
		$entrez_ids{$pathway}->[$i] = $ids[$i];
	}
}
close(E);

# read input data, write output
while(<>)
{
	chomp;
	my @fields = split /\t/;
	my $pathway = shift(@fields);
	$pathway =~ s/\)//g;	
	my $url = shift(@fields);
	my @symbols = @fields;
	
	# MuSiC required format:
	# This is a tab-delimited file prepared from a pathway database (such as KEGG), with the columns: 
	# [path_id, path_name, class, gene_line, diseases, drugs, description] The latter three columns 
	# are optional (but are available on KEGG). The gene_line contains the "entrez_id:gene_name" of 
	# all genes involved in this pathway, each separated by a "|" symbol.
	# Example:
	# hsa00061      Fatty acid biosynthesis	      Lipid Metabolism     31:ACACA|32:ACACB|27349:MCAT|2194:FASN|54995:OXSM|55301:OLAH	
	print "$pathway\t";
	print "$pathway\t";
	print "$category\t";
	
	for (my $i = 0; $i < @symbols; $i++)
	{
		print $entrez_ids{$pathway}->[$i], ":";
		print $symbols[$i];
		print "|" if ($i < @symbols-1);
	}
	
	print "\t\t\t\n";
}
