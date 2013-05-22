use warnings;
use strict;

use Vcf;
use Data::Dumper;
use Getopt::Long;

my ($vcf_out, $header);
GetOptions
(
	"vcf-out=s" => \$vcf_out,  # filtered VCF output file
	"header" => \$header  # if set, write header line to output
);

my $debug = 1;

my $patient = $ARGV[0] or die "ERROR: patient not specified\n";
my $cmp_type = $ARGV[1] or die "ERROR: Comparison type ('rem_dia' or 'rem_rel') not specified\n";
die "ERROR: invalid comparison type: $cmp_type\n" if ($cmp_type ne 'rem_dia' and $cmp_type ne 'rem_rel');
my $vcf_file = $ARGV[2] or die "ERROR: VCF file not specified\n";
my $rem_sample = $ARGV[3] or die "ERROR: remission sample name not specified\n";
my $cmp_sample = $ARGV[4] or die "ERROR: comparator sample name (diagnosis or remission) not specified\n";
my $var_type = $ARGV[5] or die "ERROR: variant type not specified ('snp' or 'indel')\n";
die "ERROR: invalid variant type: $var_type\n" if ($var_type ne 'snp' and $var_type ne 'indel');

$| = 1; # turn on autoflush

my %variant_stati = 
(
	0 => 'wildtype',
	1 => 'germline',
	2 => 'somatic',
	3 => 'LOH',
	4 => 'post-transcriptional modification',
	5 => 'unknown'
);

print STDERR "Processing file $vcf_file...\n";

my $vcf = Vcf->new(file => "$vcf_file");
$vcf->parse_header();
my (@samples) = $vcf->get_samples();

if ($vcf_out) 
{
	my $cmd = "grep -P '^#' $vcf_file > $vcf_out";
	system($cmd) == 0 or die "ERROR: grep vcf header failed: $cmd\n";
	open(VCFOUT,">>$vcf_out") or die "ERROR: could not write to file $vcf_out\n";
}

# sanity check
die "ERROR: Sample name $rem_sample not found!\n" if ($rem_sample ne $samples[0] and $rem_sample ne $samples[1]);
die "ERROR: Sample name $cmp_sample not found!\n" if ($cmp_sample ne $samples[0] and $cmp_sample ne $samples[1]);

my ($tot_var, $filtered_qual, $filtered_gt, $filtered_alt, $filtered_germ) = (0, 0, 0, 0, 0);
my %qual_num;

if ($header)
{
	print "patient\t";		
	print "sample\t";
	print "var_type\t";
	print "chr\t";
	print "pos\t";
	print "dbSNP\t";
	print "ref\t";
	print "alt\t";
	print "gene\t";
	print "add_genes\t";
	print "impact\t";
	print "effect\t";
	#print "variant_status\t";
	print "depth_rem\t";
	print "depth_leu\t";
	print "freq\t";
	print "effect\t";
	print "\n";	
}
	
while (my $line = $vcf->next_line())
{
	my $x = $vcf->next_data_hash($line);

	$tot_var ++;
	$qual_num{$x->{FILTER}->[0]} = $qual_num{$x->{FILTER}->[0]} ? $qual_num{$x->{FILTER}->[0]} + 1 : 1;
	
	if ($x->{gtypes}{$rem_sample}{GT} eq $x->{gtypes}{$cmp_sample}{GT}) # no difference in genotype?
	{
		$filtered_gt ++;
		next;
	}
			
	my $gt_rem = $x->{gtypes}{$rem_sample}{GT};
	die "ERROR: Could not determine genotype of sample $rem_sample in file $vcf_file\n" if (!defined $gt_rem or $gt_rem eq "");

	if ($gt_rem =~ /1/) # germline variant?
	{
		$filtered_germ ++;
		next;
	}
	if (@{$x->{ALT}} != 1) # more than one alternative allele?
	{
		$filtered_alt ++;
		next;
	}		

	if ($x->{FILTER}->[0] eq "REJECT") # quality filter
	{
		$filtered_qual ++;
		next;
	}
	
	if ($var_type eq 'indel') # indel
	{
		# insufficient read depth
		if ($x->{INFO}{T_DP} < 10)
		{
			print STDERR "REJECT: READ DEPTH < 10: ", $x->{INFO}{T_DP},"\n";
			next;			
		}
		
		# require support from both strands
		##INFO=<ID=T_SC,Number=4,Type=Integer,Description="In TUMOR: strandness: counts of forward-/reverse-aligned indel-supporting reads / forward-/reverse-aligned reference supporting reads">
		my ($reads_indel_fwd, $reads_indel_rev, $reads_ref_fwd, $reads_ref_rev) = split(",", $x->{INFO}{T_SC});
		if ($reads_indel_fwd == 0 or $reads_indel_rev == 0)
		{
			print STDERR "REJECT: STRAND BIAS: ",$x->{CHROM},":",$x->{POS},"\t","reads_indel_fwd: $reads_indel_fwd\treads_indel_rev: $reads_indel_rev\n";
			next;
		}
		
		# require high consensus call for indel
		##INFO=<ID=T_AC,Number=2,Type=Integer,Description="In TUMOR: # of reads supporting consensus indel/any indel at the site">
		my ($num_reads_cons_indel, $num_reads_any_indel) = split(",", $x->{INFO}{T_AC});
		if ($num_reads_cons_indel/$num_reads_any_indel < 0.7)
		{
			print STDERR "REJECT: BAD CONSENSUS: ",$x->{CHROM},":",$x->{POS},"\n";
			next;
		}
		
		# check alignment quality around indel
		##INFO=<ID=T_NQSMM,Number=2,Type=Float,Description="In TUMOR: Within NQS window: fraction of mismatching bases in consensus indel-supporting reads/in reference-supporting reads">
		my ($frac_mm_reads_indel, $frac_mm_reads_ref) = split(",", $x->{INFO}{T_NQSMM});
		if ($frac_mm_reads_indel > $frac_mm_reads_ref)
		{
			print STDERR "REJECT: POOR ALIGNMENT: ",$x->{CHROM},":",$x->{POS},"\t","frac_mm_reads_indel: $frac_mm_reads_indel\tfrac_mm_reads_ref: $frac_mm_reads_ref\n";
			next;
		}
		
				
#		print "reads_indel_fwd: $reads_indel_fwd\n";
#		print "reads_indel_rev: $reads_indel_rev\n";
#		print "reads_ref_fwd: $reads_ref_fwd\n";
#		print "reads_ref_rev: $reads_ref_rev\n";
	}

	print VCFOUT "$line" if ($vcf_out);
	
	print "$patient\t";		
	print "$cmp_type\t";
	print "$var_type\t";
	print $x->{CHROM},"\t";
	print $x->{POS},"\t";
	print $x->{ID},"\t";
	print $x->{REF},"\t";
	print $x->{ALT}->[0],"\t";
	my ($gene, $add_genes, $impact, $effect) = get_impact($x->{INFO}{EFF});
	print "$gene\t";
	print "$add_genes\t";
	print "$impact\t";
	print "$effect\t";
#	print join(",", @{$x->{FILTER}}),"\t";
#	print exists $x->{gtypes}{$cmp_sample}{SS} ? $variant_stati{$x->{gtypes}{$cmp_sample}{SS}} : "n/a", "\t";
	if ($var_type eq 'snp')
	{
		print exists $x->{gtypes}{$rem_sample}{DP} ? $x->{gtypes}{$rem_sample}{DP} : "n/a", "\t";
		print exists $x->{gtypes}{$cmp_sample}{DP} ? $x->{gtypes}{$cmp_sample}{DP} : "n/a", "\t";
		print exists $x->{gtypes}{$cmp_sample}{FA} ? $x->{gtypes}{$cmp_sample}{FA} : "n/a", "\t";	
	}
	else
	{
		my ($depth_rem, $depth_leu) = ($x->{INFO}{N_DP}, $x->{INFO}{T_DP});
		print defined $depth_rem ? $depth_rem : "n/a", "\t";
		print defined $depth_leu ? $depth_leu : "n/a", "\t";
		
		##INFO=<ID=T_AC,Number=2,Type=Integer,Description="In TUMOR: # of reads supporting consensus indel/any indel at the site">
		my ($num_reads_cons_indel, $num_reads_any_indel) = split(",", $x->{INFO}{T_AC});
		my $all_frac = $depth_leu ? $num_reads_cons_indel / $depth_leu : "n/d";
		print "$all_frac\t";
	}
	print "EFF=",$x->{INFO}{EFF},"\t";
	print "\n";
		
#	print "\n"; print Dumper($x); exit;
}
$vcf->close();
close(VCFOUT) if ($vcf_out);
	
if ($debug)
{
	print STDERR "  Total number of variants: $tot_var\n";
	print STDERR "  Variants by quality:\n";
	foreach my $k (keys(%qual_num))
	{
		print STDERR "    $k: ", $qual_num{$k},"\n";
	}
	print STDERR "  Excluded due to quality: $filtered_qual\n";
	print STDERR "  Excluded due to equal genotype: $filtered_gt\n";
	print STDERR "  Excluded due to missing alternative allele: $filtered_alt\n";
	print STDERR "  Excluded germline variants: $filtered_germ\n";
}

# ------------------------------------------

sub get_impact
{
	my $effs = shift or die "ERROR: effect not specified";

	# determine all genes impacted by variants
	my (%genes_by_impact, %all_genes, $combined_effect, $combined_impact);	
	foreach my $eff (split(",", $effs))
	{
		my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
			or die "ERROR: could not parse SNP effect: $effs\n";

		my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
			$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
				or die "ERROR: could not parse SNP effect: $eff\n"; 

		# gene impacted by variant?
		if ($gene_name)
		{
			$genes_by_impact{$impact}{$gene_name} = $effect;
			$all_genes{$gene_name} = 1;
		}

		$combined_impact = $impact;		
		$combined_effect = $effect;
	}
	
	# if multiple genes are affected, preferentially chose gene with the predicted higher impact
	if ($genes_by_impact{HIGH})
	{
		$combined_impact = "HIGH";
	}
	elsif ($genes_by_impact{MODERATE})
	{
		$combined_impact = "MODERATE";
	}
	elsif ($genes_by_impact{LOW})
	{
		$combined_impact = "LOW";
	}
	elsif ($genes_by_impact{MODIFIER})
	{
		$combined_impact = "MODIFIER";
	}
	
	my ($gene, $add_genes) = ("", "");
	if (keys(%all_genes) > 0)
	{
		my @sorted_genes = sort keys(%{$genes_by_impact{$combined_impact}});
		$gene = $sorted_genes[0]; # first choice is first in alphabetically sorted list
		if ($gene =~ /^LOC/) # if this is a generic gene name, try to find non-generic one instead
		{
			foreach my $g (@sorted_genes)
			{
				if ($g !~ /^LOC/)
				{
					$gene = $g;
					last;
				}	
			}
		}
		$combined_effect = $genes_by_impact{$combined_impact}{$gene};
		delete $all_genes{$gene};
		$add_genes = join(",", keys(%all_genes));
	}
#		# determine overall impact
#		if ($combined_impact eq "n/d" or $combined_impact eq "MODIFIER")
#		{
#			$gene = $gene_name;
#			$combined_impact = $impact;
#			$combined_effect = $effect;
#		}
#		elsif ($combined_impact eq "LOW" and $impact =~ /(MODERATE|HIGH)/)
#		{
#			$gene = $gene_name;
#			$combined_impact = $impact;
#			$combined_effect = $effect;
#		}
#		elsif ($combined_impact eq "MODERATE" and $impact =~ /HIGH/)
#		{
#			$gene = $gene_name;
#			$combined_impact = $impact;
#			$combined_effect = $effect;
#		}
#		elsif ($impact eq "HIGH")
#		{
#			$gene = $gene_name;
#			$combined_impact = $impact;
#			$combined_effect = $effect;
#		}		

	return ($gene, $add_genes, $combined_impact, $combined_effect);
}
