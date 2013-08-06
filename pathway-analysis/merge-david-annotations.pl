# takes results from the two david runs (using different input order of genes) and merges them into a final result
# necessarry because david web interface appears to not output annotations for the first gene of the list

use warnings FATAL => qw( all );
use strict;

use Carp;

my $file1 = $ARGV[0] or die "file 1 not specified\n";
my $file2 = $ARGV[1] or die "file 2 not specified\n";

my %p2g;
open(S, $file1) or die "error opening file\n";
while(<S>)
{
	chomp;
	my ($id, $name, $class, $gene_line) = split(/\t/);
	map { $p2g{$id}{$_} = 1 } split(/\|/, $gene_line);
}
close(S);

my %added_genes;
my $added_annotations = 0;
my $single = 0;
open(S, $file2);
while(<S>)
{
	my ($id, $name, $class, $gene_line) = split(/\t/) or die "error opening file\n";
	
	my @genes = keys(%{$p2g{$id}});
	map 
	{ 
		if (!exists $p2g{$id}{$_}) 
		{ 
#			print "$id: $_\n"; 
			push(@genes, $_); 
			$added_annotations ++;
			$added_genes{$_} = 1;
		} 
	} split(/\|/, $gene_line);

#	if (@genes > 1000)
#	{
#		print STDERR "Excluding pathway $id:$name because of too many genes (".scalar(@genes).")";
#		next;
#	}
	next if ($name eq "ultra-high-sulfur keratin");
	next if ($name eq "complete proteome");
#	next if ($name eq "ribonucleotide binding");
#	next if ($name eq "regulation of RNA metabolic process");
#	next if ($name eq "regulation of transcription, DNA-dependent");
#	next if ($name eq "metal ion binding");
	next if ($class eq "SP_PIR_KEYWORDS");
	next if ($class eq "UP_SEQ_FEATURE");
	next if ($class eq "GOTERM_CC_FAT");
	
#	if (@genes < 2)
#	{
#		$single ++;
#		next;
#	}	

	print "$id\t$name\t$class\t",join("\|", @genes),"\t\t\t\n";
}
close(S);

print STDERR "Added ",scalar(keys(%added_genes))," genes: ",join(",", keys(%added_genes)),"\n";
print STDERR "Addded $added_annotations annotations.\n";
print STDERR "Excluded $single pathways because they contain only a single gene.\n" if ($single > 0);