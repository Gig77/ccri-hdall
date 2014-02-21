use warnings FATAL => qw( all );
use strict;

use Carp;
use Getopt::Long;

# normalize allelic frequencies based on copy-number

my ($cnv_file);
GetOptions
(
	"cnv-file=s" => \$cnv_file
);

croak "ERROR: --cnv-file not specified" if (!$cnv_file);

# TABLE: hdall.cnv.tsv
my %cnv;
my $read_lines = 0;
open(CNV,"$cnv_file") or croak "could not open file $cnv_file";
<CNV>; # skip header
while(<CNV>)
{
	chomp;
	my ($patient, $sample, $event, $cnumber, $chromosome, $start, $end, $tsize, $num_genes, $genes) = split(/\t/);
	$chromosome = "chr$chromosome";
	
	$cnv{"$patient\t$sample"} = 1;
	$cnv{"$patient\t$sample\t$chromosome"}{"$start-$end"} = $cnumber;
	$read_lines ++;
}
close(CNV);
print STDERR "INFO: $read_lines lines read from file $cnv_file\n";

# TABLE: filtered-variants.cosmic
my $normalized = 0;
my $header = <>;
chomp($header);
print "$header\tcopy_no\tfreq_leu_norm\n";
while(<>)
{
	chomp;
	my ($patient, $sample, $var_type, $status, $rejected_because, $chr, $pos, $dbSNP, $ref, $alt, $gene, $add_genes, $impact, $effect, $non_silent, $deleterious, $exons, 
		$dp_rem_tot, $dp_rem_ref, $dp_rem_var, $freq_rem, $dp_leu_tot, $dp_leu_ref, $dp_leu_var, $freq_leu) = split("\t");
	my $line = $_;
	
	$sample =~ s/rem_rel3/rem_rel/; # normalize relapse 3 of patient 715

	my ($copy_no, $norm_af) = ("n/a", $freq_leu);
	if ($cnv{"$patient\t$sample"})
	{
		$copy_no = 2;  # default
		
		if ($cnv{"$patient\t$sample\t$chr"})
		{
			foreach my $coord (keys(%{$cnv{"$patient\t$sample\t$chr"}}))
			{
				my ($start, $end) = split("-", $coord);
				if ($pos >= $start and $pos <= $end)
				{
					$copy_no = $cnv{"$patient\t$sample\t$chr"}{$coord};
					if ($cnv{"$patient\t$sample\t$chr"}{$coord} > 2)
					{
						$norm_af = $freq_leu / 2 * $copy_no;
						$normalized ++;			
					}
					last;
				}
			}		
		}		
	}
	
	print join("\t", $line), "\t$copy_no\t$norm_af\n";
}
print STDERR "INFO: Allelic frequency of $normalized variants adjusted to copy number > 2.\n";
