use warnings FATAL => qw( all );
use strict;

use Carp;
use Getopt::Long;

# STDIN: CNV data (hdall.cnv.tsv)
# STDOUT: CNVs in MAF format such that it can be used with music pathscan

my ($max_size_bp);
GetOptions
(
	"max-size-bp=i" => \$max_size_bp
);

croak ("ERROR: --max-size-bp not specified\n") if (!defined $max_size_bp);

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
print STDERR "INFO: ", (scalar(keys(%id2sym)), " id mappgins read from file $ENV{HOME}/hdall/results/id-mappings.tsv\n");


# OUTPUT format
#    genome music path-scan --help:
#    	Only the following four columns in the MAF are used. All other columns may be left blank.
#
#     	Col 1: Hugo_Symbol (Need not be HUGO, but must match gene names used in the pathway file)
#     	Col 2: Entrez_Gene_Id (Matching Entrez ID trump gene name matches between pathway file and MAF)
#     	Col 9: Variant_Classification
#     	Col 16: Tumor_Sample_Barcode (Must match the name in sample-list, or contain it as a substring)
#
#   	The Entrez_Gene_Id can also be left blank (or set to 0), but it is highly recommended, in case
#	    genes are named differently in the pathway file and the MAF file.

# TABLE: cnv.maf
<>; # skip header
while(<>)
{
	chomp;
	my ($patient, $sample, $event, $cnumber, $chromosome, $start, $end, $size, $num_genes, $genestr) = split("\t");
	next if ($patient =~ /^(1187|27|389|417|53|54|721|826)$/);
	next if ($sample =~ /rem_rel3/); # otherwise gives PathScan error b/c of missing WIG file
	$sample =~ s/rem_//;
	$genestr =~ s/\"//;
	
	if ($size <= $max_size_bp)
	{
		foreach my $gene (split(", ", $genestr))
		{
			next if (!$gene or $gene eq " ");
			$gene = $id2sym{$gene} if ($id2sym{$gene});  # remap symbol
			 
			print "$gene\t0\tArray\t";
			print "\t\t\t\t\t";
			print "".($event eq "gain" ? "Frame_Shift_Ins" : "Frame_Shift_Del")."\t";
			print "".($event eq "gain" ? "ins" : "del")."\t";
			print"\t\t\t\t\t";
			print "$patient"."_$sample\n";
		}		
	}
}

