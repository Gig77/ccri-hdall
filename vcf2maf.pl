use warnings;
use strict;

use Vcf;
use Tabix;
use Getopt::Long;
use Data::Dumper;

my $music_roi;
GetOptions
(
	"music-roi=s" => \$music_roi  # MuSiC region of interest (ROI) file (must be tabix accessible, i.e. compressed and indexed)
);

# NOTE: use filtered VCF to run this script! 
# write filtered VCF fle using the filter-variants.pl script!

die "ERROR: --music-roi not specified (.gz file)\n" if (!$music_roi);

my $roi = Tabix->new(-data => "$music_roi");

my $gene = roi2gene("chr1", 701707, 701708);
print "$gene\n";

my $vcf = Vcf->new(file => "-");

$vcf->parse_header();
my (@samples) = $vcf->get_samples();

print_header();

while (my $x = $vcf->next_data_hash())
{
	my $chr = $x->{CHROM};
	my $pos = $x->{POS};
	my $gene = roi2gene($chr, $pos, $pos+1);

	print "$chr\t$pos\t$gene\n";		 

#	my ($genes, $region, $impact) = get_impact($x->{INFO}{EFF});
}
$vcf->close();


#while (my $x = $vcf->next_data_hash())
#{	 
#	$tot_var ++;
#	$qual_num{$x->{FILTER}->[0]} = $qual_num{$x->{FILTER}->[0]} ? $qual_num{$x->{FILTER}->[0]} + 1 : 1;
#	
#	if ($x->{gtypes}{$rem_sample}{GT} eq $x->{gtypes}{$cmp_sample}{GT}) # no difference in genotype?
#	{
#		$filtered_gt ++;
#		next;
#	}
#			
#	my $gt_rem = $x->{gtypes}{$rem_sample}{GT};
#	die "ERROR: Could not determine genotype of sample $rem_sample in file $vcf_file\n" if (!defined $gt_rem or $gt_rem eq "");
#
#	if ($gt_rem =~ /1/) # germline variant?
#	{
#		$filtered_germ ++;
#		next;
#	}
#	if (@{$x->{ALT}} != 1) # more than one alternative allele?
#	{
#		$filtered_alt ++;
#		next;
#	}		
#
#	if ($x->{FILTER}->[0] eq "REJECT") # quality filter
#	{
#		$filtered_qual ++;
#		next;
#	}
#	
#	if ($var_type eq 'indel') # indel
#	{
#		# insufficient read depth
#		if ($x->{INFO}{T_DP} < 10)
#		{
#			print STDERR "REJECT: READ DEPTH < 10: ", $x->{INFO}{T_DP},"\n";
#			next;			
#		}
#		
#		# require support from both strands
#		##INFO=<ID=T_SC,Number=4,Type=Integer,Description="In TUMOR: strandness: counts of forward-/reverse-aligned indel-supporting reads / forward-/reverse-aligned reference supporting reads">
#		my ($reads_indel_fwd, $reads_indel_rev, $reads_ref_fwd, $reads_ref_rev) = split(",", $x->{INFO}{T_SC});
#		if ($reads_indel_fwd == 0 or $reads_indel_rev == 0)
#		{
#			print STDERR "REJECT: STRAND BIAS: ",$x->{CHROM},":",$x->{POS},"\t","reads_indel_fwd: $reads_indel_fwd\treads_indel_rev: $reads_indel_rev\n";
#			next;
#		}
#		
#		# require high consensus call for indel
#		##INFO=<ID=T_AC,Number=2,Type=Integer,Description="In TUMOR: # of reads supporting consensus indel/any indel at the site">
#		my ($num_reads_cons_indel, $num_reads_any_indel) = split(",", $x->{INFO}{T_AC});
#		if ($num_reads_cons_indel/$num_reads_any_indel < 0.7)
#		{
#			print STDERR "REJECT: BAD CONSENSUS: ",$x->{CHROM},":",$x->{POS},"\n";
#			next;
#		}
#		
#		# check alignment quality around indel
#		##INFO=<ID=T_NQSMM,Number=2,Type=Float,Description="In TUMOR: Within NQS window: fraction of mismatching bases in consensus indel-supporting reads/in reference-supporting reads">
#		my ($frac_mm_reads_indel, $frac_mm_reads_ref) = split(",", $x->{INFO}{T_NQSMM});
#		if ($frac_mm_reads_indel > $frac_mm_reads_ref)
#		{
#			print STDERR "REJECT: POOR ALIGNMENT: ",$x->{CHROM},":",$x->{POS},"\t","frac_mm_reads_indel: $frac_mm_reads_indel\tfrac_mm_reads_ref: $frac_mm_reads_ref\n";
#			next;
#		}
#		
#				
##		print "reads_indel_fwd: $reads_indel_fwd\n";
##		print "reads_indel_rev: $reads_indel_rev\n";
##		print "reads_ref_fwd: $reads_ref_fwd\n";
##		print "reads_ref_rev: $reads_ref_rev\n";
#	}
#
#	print "$patient\t";		
#	print "$cmp_type\t";
#	print "$var_type\t";
#	print $x->{CHROM},"\t";
#	print $x->{POS},"\t";
#	print $x->{ID},"\t";
#	print $x->{REF},"\t";
#	print $x->{ALT}->[0],"\t";
#	my ($genes, $region, $impact) = get_impact($x->{INFO}{EFF});
#	print "$genes\t";
#	print "$impact\t";
##	print "$region\t";
##	print join(",", @{$x->{FILTER}}),"\t";
##	print exists $x->{gtypes}{$cmp_sample}{SS} ? $variant_stati{$x->{gtypes}{$cmp_sample}{SS}} : "n/a", "\t";
#	if ($var_type eq 'snp')
#	{
#		print exists $x->{gtypes}{$rem_sample}{DP} ? $x->{gtypes}{$rem_sample}{DP} : "n/a", "\t";
#		print exists $x->{gtypes}{$cmp_sample}{DP} ? $x->{gtypes}{$cmp_sample}{DP} : "n/a", "\t";
#		print exists $x->{gtypes}{$cmp_sample}{FA} ? $x->{gtypes}{$cmp_sample}{FA} : "n/a", "\t";	
#	}
#	else
#	{
#		my ($depth_rem, $depth_leu) = ($x->{INFO}{N_DP}, $x->{INFO}{T_DP});
#		print defined $depth_rem ? $depth_rem : "n/a", "\t";
#		print defined $depth_leu ? $depth_leu : "n/a", "\t";
#		
#		##INFO=<ID=T_AC,Number=2,Type=Integer,Description="In TUMOR: # of reads supporting consensus indel/any indel at the site">
#		my ($num_reads_cons_indel, $num_reads_any_indel) = split(",", $x->{INFO}{T_AC});
#		my $all_frac = $depth_leu ? $num_reads_cons_indel / $depth_leu : "n/d";
#		print "$all_frac\t";
#	}
#	print "EFF=",$x->{INFO}{EFF},"\t";
#	print "\n";
#		
##	print "\n"; print Dumper($x); exit;
#}
#$vcf->close();
#	
#if ($debug)
#{
#	print STDERR "  Total number of variants: $tot_var\n";
#	print STDERR "  Variants by quality:\n";
#	foreach my $k (keys(%qual_num))
#	{
#		print STDERR "    $k: ", $qual_num{$k},"\n";
#	}
#	print STDERR "  Excluded due to quality: $filtered_qual\n";
#	print STDERR "  Excluded due to equal genotype: $filtered_gt\n";
#	print STDERR "  Excluded due to missing alternative allele: $filtered_alt\n";
#	print STDERR "  Excluded germline variants: $filtered_germ\n";
#}
#
## ------------------------------------------
#

sub print_header
{
	print "#version 2.3\n";
	
	print "Hugo_Symbol\t";
	print "Entrez_Gene_Id\t";
	print "Center\t";
	print "NCBI_Build\t";
	print "Chromosome\t";
	print "Start_Position\t";
	print "End_Position\t";
	print "Strand\t";
	print "Variant_Classification\t";
	print "Variant_Type\t";
	print "Reference_Allele\t";
	print "Tumor_Seq_Allele1\t";
	print "Tumor_Seq_Allele2\t";
	print "dbSNP_RS\t";
	print "dbSNP_Val_Status\t";
	print "Tumor_Sample_Barcode\t";
	print "Matched_Norm_Sample_Barcode\t";
	print "Match_Norm_Seq_Allele1\t";
	print "Match_Norm_Seq_Allele2\t";
	print "Tumor_Validation_Allele1\t";
	print "Tumor_Validation_Allele2\t";
	print "Match_Norm_Validation_Allele1\t";
	print "Match_Norm_Validation_Allele2\t";
	print "Verification_Status\t";
	print "Validation_Status\t";
	print "Mutation_Status\t";
	print "Sequencing_Phase\t";
	print "Sequence_Source\t";
	print "Validation_Method\t";
	print "Score\t";
	print "BAM_File\t";
	print "Sequencer\t";
	print "Tumor_Sample_UUID\t";
	print "Matched_Norm_Sample_UUID\n";	
}

sub get_impact
{
	my $effs = shift or die "ERROR: effect not specified";

	my %genes;
	my $overall_impact = "n/d";
	foreach my $eff (split(",", $effs))
	{
		my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
			or die "ERROR: could not parse SNP effect: $effs\n";

		my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
			$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
				or die "ERROR: could not parse SNP effect: $eff\n"; 

		# determine overall impact
		if ($overall_impact eq "n/d" or $overall_impact eq "MODIFIER")
		{
			$overall_impact = $impact;
		}
		elsif ($overall_impact eq "LOW" and $impact =~ /(MODERATE|HIGH)/)
		{
			$overall_impact = $impact;
		}
		elsif ($overall_impact eq "MODERATE" and $impact =~ /HIGH/)
		{
			$overall_impact = $impact;
		}
		elsif ($impact eq "HIGH")
		{
			$overall_impact = $impact;
		}

		$genes{$gene_name} = $effect
			if ($gene_name);
	}
	
#	my $impact = "";
#	my $region = "";
#	foreach my $g (keys(%genes))
#	{
#		$impact .= $g;
#		$region .= $genes{$g}.")";
#	}
	
	return (join(",", keys(%genes)), join(",", values(%genes)), $overall_impact);
}

#sub roi2gene
#{
#	my ($chr, $start, $end) = @_;
#	my $cmd = "echo $'$chr\t$start\t$end' | /home/STANNA/christian.frech/tools/bedtools-2.17.0/bin/intersectBed -a ~/hdall/results/music/ucsc-genes.hg19.roi -b stdin 
#	my $result = ` 
#	
#}

sub roi2gene
{
	my ($chr, $start, $end) = @_;
	my $iter = $roi->query("chr1", $start, $end);
	my $line = $roi->read($iter);
	
	if ($line)
	{
		return (split("\t", $line))[3]; 	
	}
	else
	{
		return "";
	}
}