use warnings FATAL => qw( all );
use strict;

use lib "$ENV{HOME}/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Vcf;
use Data::Dumper;
use Getopt::Long;
use Tabix;
use Carp;

my ($vcf_out, $header, $rejected_variants_file, $sample_identifier, $vcf_in, $var_type);
my ($rmsk_file, $simplerepeat_file, $blacklist_file, $segdup_file, $g1k_accessible_file, $ucsc_retro_file, $remission_variants_file);
GetOptions
(
	"sample=s" => \$sample_identifier, # e.g. 314_rem_dia
	"vcf-in=s" => \$vcf_in,  # VCF input file
	"vcf-out=s" => \$vcf_out,  # filtered VCF output file
	"variant-type=s" => \$var_type,  # 'snp' or 'indel'
	"header" => \$header,  # if set, write header line to output
	"rmsk-file=s" => \$rmsk_file, # TABIX indexed UCSC table rmsk
	"simpleRepeat-file=s" => \$simplerepeat_file, # TABIX indexed UCSC table rmsk
	"blacklist-file=s" => \$blacklist_file, # TABIX indexed UCSC table wgEncodeDacMapabilityConsensusExcludable
	"segdup-file=s" => \$segdup_file, # TABIX indexed UCSC table genomicSuperDups
	"g1k-accessible=s" => \$g1k_accessible_file, # TABIX indexed UCSC table tgpPhase1AccessibilityPilotCriteria
	"ucscRetro=s" => \$ucsc_retro_file, # TABIX indexed UCSC table ucscRetroAli5
	"rejected-variants-file=s" => \$rejected_variants_file, # file with variants rejected based on manual curation; will be filtered from output
	"remission-variants-file=s" => \$remission_variants_file # TABIX indexed file with variants found in remission samples (GATK)
);

# TABLE: filtered-variants
if ($header)
{
	print "patient\t";		
	print "sample\t";
	print "var_type\t";
	print "status\t";
	print "rejected_because\t";
	print "chr\t";
	print "pos\t";
	print "dbSNP\t";
	print "ref\t";
	print "alt\t";
	print "gene\t";
	print "add_genes\t";
	print "impact\t";
	print "effect\t";
	print "exons\t";
	print "dp_rem_tot\t";
	print "dp_rem_ref\t";
	print "dp_rem_var\t";
	print "freq_rem\t";
	print "dp_leu_tot\t";
	print "dp_leu_ref\t";
	print "dp_leu_var\t";
	print "freq_leu\t";
	print "aa_change\t";
	print "effect\t";
	print "Polyphen2\t";
	print "SIFT\t";
	print "GERP++\t";
	print "SiPhy\t";
	print "InterPro\t";
	print "AF_1000G\t";
	print "repeat\t";
	print "segdup\t";
	print "blacklist\t";
	print "g1k-accessible\t";	
	print "retro\t";
	print "rem_samples\n";	
	exit;
}

my $debug = 1;

croak "ERROR: --sample not specified" if (!$sample_identifier);
croak "ERROR: --vcf-in not specified" if (!$vcf_in);
croak "ERROR: --variant-type not specified ('snp' or 'indel')" if (!$var_type);
croak "ERROR: invalid variant type: $var_type\n" if ($var_type ne 'snp' and $var_type ne 'indel');
croak "ERROR: --rmsk-file not specified" if (!$rmsk_file);
croak "ERROR: --simpleRepeat-file not specified" if (!$simplerepeat_file);
croak "ERROR: --blacklist-file not specified" if (!$blacklist_file);
croak "ERROR: --segdup-file not specified" if (!$segdup_file);
croak "ERROR: --g1k-accessible not specified" if (!$g1k_accessible_file);
croak "ERROR: --ucscRetro not specified" if (!$ucsc_retro_file);
croak "ERROR: --remission-variants-file not specified" if (!$remission_variants_file);

my ($patient, $rem_sample, $cmp_sample) = split("_", $sample_identifier) or croak "ERROR: could not parse sample identifier\n";
my $cmp_type = $rem_sample."_".$cmp_sample;
die "ERROR: invalid comparison type: $cmp_type\n" if ($cmp_type ne 'rem_dia' and $cmp_type ne 'rem_rel');

my %patient2sample = (
	'A_rem' => 'A13324_rem',
	'A_dia' => 'A12642_dia',
	'A_rel' => 'A12886_rel',
	'B_rem' => 'B20946_rem',
	'B_dia' => 'B19668_dia',
	'B_rel' => 'B15010_rel',
	'C_rem' => 'C20499_rem',
	'C_dia' => 'C19797_dia',
	'C_rel' => 'C15050_rel',
	'D_rem' => 'D4502_rem',
	'D_dia' => 'D3826_dia',
	'D_rel' => 'D10183_rel',
	'E_rem' => 'E13861_rem',
	'E_dia' => 'E13174_dia',
	'E_rel' => 'E13479_rel',
	'X_rem' => 'X1847_rem',
	'X_dia' => 'X1286_dia',
	'X_rel' => 'X12831_rel',
	'Y_rem' => 'Y3767_rem',
	'Y_dia' => 'Y3141_dia',
	'Y_rel' => 'Y10284_rel'
);

$rem_sample = $patient2sample{$patient."_$rem_sample"} ? $patient2sample{$patient."_$rem_sample"} : $patient."_$rem_sample"; 
$cmp_sample = $patient2sample{$patient."_$cmp_sample"} ? $patient2sample{$patient."_$cmp_sample"} : $patient."_$cmp_sample"; 

# read kgXref, knownCanonical to determine UCSC canonical transcripts affected by variant
my %kgID2refSeq;
open(G,"$ENV{HOME}/generic/data/hg19/hg19.kgXref.txt") or die "could not open file $ENV{HOME}/generic/data/hg19/hg19.kgXref.txt";
while(<G>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refSeq, $protAcc, $description, $rfamAcc, $tRnaName) = split(/\t/);

	$kgID2refSeq{$kgID} = $refSeq if ($refSeq);
}
close(G);
INFO(scalar(keys(%kgID2refSeq))." gene descriptions read from file $ENV{HOME}/generic/data/hg19/hg19.kgXref.txt");

my %canonical;
open(G,"$ENV{HOME}/generic/data/hg19/hg19.knownCanonical.txt") or die "could not open file $ENV{HOME}/generic/data/hg19/hg19.knownCanonical.txt";
<G>; # skip header
while(<G>)
{
	chomp;
	my ($chrom, $chromStart, $chromEnd, $clusterId, $transcript, $protein) = split(/\t/);
	
	$canonical{$kgID2refSeq{$transcript}} = 1 if ($kgID2refSeq{$transcript});
}
close(G);
INFO(scalar(keys(%canonical))." canonical genes read from file $ENV{HOME}/generic/data/hg19/hg19.knownCanonical.txt");

my %rejected_variants;
if ($rejected_variants_file)
{
	open(R,"$rejected_variants_file") or die "could not open file $rejected_variants_file";
	<R>; # skip header
	while(<R>)
	{
		chomp;
		my ($patient, $sample, $var_type, $rejected_because, $chr, $pos) = split(/\t/);
		$rejected_variants{"$patient\t$sample\t$chr\t$pos"} = $rejected_because;
	}
	close(R);
	INFO(scalar(keys(%rejected_variants))." variants read from file $rejected_variants_file");
}

my $rmsk = Tabix->new(-data => $rmsk_file);
my $simpleRepeat = Tabix->new(-data => $simplerepeat_file);
my $blacklistdb = Tabix->new(-data => $blacklist_file);
my $segdup = Tabix->new(-data => $segdup_file);
my $g1kAcc = Tabix->new(-data => $g1k_accessible_file);
my $ucscRetro = Tabix->new(-data => $ucsc_retro_file);
my $remission = Tabix->new(-data => $remission_variants_file);

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

INFO("Processing file $vcf_in...");

my $vcf = Vcf->new(file => "$vcf_in");
$vcf->parse_header();

my $mutect = $vcf->get_header_line(key => 'GATKCommandLine', ID => 'MuTect')->[0]->{'CommandLineOptions'};
$mutect = $vcf->get_header_line(key => 'MuTect')->[0]->[0]->{'value'} if (!$mutect);
if ($mutect)
{
	($rem_sample) = $mutect =~ /normal_sample_name=(\S+)/;
	($cmp_sample) = $mutect =~ /tumor_sample_name=(\S+)/;
}
print STDERR "Normal sample name: $rem_sample\n";
print STDERR "Tumor sample name: $cmp_sample\n";

my (@samples) = $vcf->get_samples();

if ($vcf_out) 
{
	my $cmd = "grep -P '^#' $vcf_in > $vcf_out";
	system($cmd) == 0 or die "ERROR: grep vcf header failed: $cmd\n";
	open(VCFOUT,">>$vcf_out") or die "ERROR: could not write to file $vcf_out\n";
}

# sanity check
die "ERROR: Sample name $rem_sample not found!\n" if ($rem_sample ne $samples[0] and $rem_sample ne $samples[1]);
die "ERROR: Sample name $cmp_sample not found!\n" if ($cmp_sample ne $samples[0] and $cmp_sample ne $samples[1]);

my ($tot_var, $filtered_qual, $filtered_gt, $filtered_alt, $filtered_germ) = (0, 0, 0, 0, 0);
my ($numrep, $num_blacklist, $numsegdup, $num_not_accessible, $num_retro, $num_remission) = (0, 0, 0, 0, 0, 0);
my %qual_num;

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
	die "ERROR: Could not determine genotype of sample $rem_sample in file $vcf_in\n" if (!defined $gt_rem or $gt_rem eq "");

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

	my $status = $x->{FILTER}->[0];
	$status = "MISSED" if ($patient eq "C" and $x->{CHROM} eq "chr16" and $x->{POS} eq "3789627"); # keep mutation CREBBP mutation falsely rejected by MuTect
	
	if ($status eq "REJECT") # rejected by MuTect
	{
		$filtered_qual ++;
		next;
	}

	my ($dp_tum, $dp_rem, $freq_tum, $freq_rem, $ad_tum_ref, $ad_tum_alt, $ad_rem_ref, $ad_rem_alt);
	
#	print Dumper($x);
#	exit;		
	
	if ($var_type eq 'snp')
	{
		##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		($ad_tum_ref, $ad_tum_alt) = split(",", $x->{gtypes}{$cmp_sample}{AD});
		($ad_rem_ref, $ad_rem_alt) = split(",", $x->{gtypes}{$rem_sample}{AD});
		
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		($dp_tum, $dp_rem) = ($x->{gtypes}{$cmp_sample}{DP}, $x->{gtypes}{$rem_sample}{DP});

		##FORMAT=<ID=FA,Number=A,Type=Float,Description="Allele fraction of the alternate allele with regard to reference">
		$freq_tum = $x->{gtypes}{$cmp_sample}{FA};
		$freq_rem = $x->{gtypes}{$rem_sample}{FA};
		
		#next if ($status eq "REJECT" and ($dp_tum <= 50 or $dp_rem <= 50 or $freq_tum < 0.2 or $freq_rem > 0.05));	
	}
	elsif ($var_type eq 'indel') # indel
	{
		##INFO=<ID=T_DP,Number=1,Type=Integer,Description="In TUMOR: total coverage at the site">
		##INFO=<ID=N_DP,Number=1,Type=Integer,Description="In NORMAL: total coverage at the site">
		($dp_tum, $dp_rem) = ($x->{INFO}{T_DP}, $x->{INFO}{N_DP});
		
		##INFO=<ID=T_AC,Number=2,Type=Integer,Description="In TUMOR: # of reads supporting consensus indel/any indel at the site">
		##INFO=<ID=N_AC,Number=2,Type=Integer,Description="In NORMAL: # of reads supporting consensus indel/any indel at the site">
		my ($ad_tum_any_indel, $ad_rem_any_indel);
		($ad_tum_alt, $ad_tum_any_indel) = split(",", $x->{INFO}{T_AC}); 		
		($ad_rem_alt, $ad_rem_any_indel) = split(",", $x->{INFO}{N_AC});
		($ad_tum_ref, $ad_rem_ref) = ($dp_tum - $ad_tum_any_indel, $dp_rem - $ad_rem_any_indel); 		
		$freq_tum = sprintf("%.3f", $ad_tum_alt / $dp_tum);
		$freq_rem = sprintf("%.3f", $ad_rem_alt / $dp_rem);

		# insufficient read depth
		if ($dp_tum < 10)
		{
			INFO("REJECT: READ DEPTH < 10: $dp_tum");
			next;			
		}

		# require high consensus call for indel
		if ($ad_tum_alt/$ad_tum_any_indel < 0.7)
		{
			INFO("REJECT: BAD CONSENSUS: ",$x->{CHROM},":",$x->{POS},"");
			next;
		}		
		
		# require support from both strands
		##INFO=<ID=T_SC,Number=4,Type=Integer,Description="In TUMOR: strandness: counts of forward-/reverse-aligned indel-supporting reads / forward-/reverse-aligned reference supporting reads">
		my ($reads_indel_fwd, $reads_indel_rev, $reads_ref_fwd, $reads_ref_rev) = split(",", $x->{INFO}{T_SC});
		if ($reads_indel_fwd == 0 or $reads_indel_rev == 0)
		{
			INFO("REJECT: STRAND BIAS: ",$x->{CHROM},":",$x->{POS},"\t","reads_indel_fwd: $reads_indel_fwd\treads_indel_rev: $reads_indel_rev");
			next;
		}
		
		# check alignment quality around indel
		##INFO=<ID=T_NQSMM,Number=2,Type=Float,Description="In TUMOR: Within NQS window: fraction of mismatching bases in consensus indel-supporting reads/in reference-supporting reads">
		my ($frac_mm_reads_indel, $frac_mm_reads_ref) = split(",", $x->{INFO}{T_NQSMM});
		if ($frac_mm_reads_indel - $frac_mm_reads_ref > 0.01)
		{
			INFO("REJECT: POOR ALIGNMENT: ",$x->{CHROM},":",$x->{POS},"\t","frac_mm_reads_indel: $frac_mm_reads_indel\tfrac_mm_reads_ref: $frac_mm_reads_ref");
			next;
		}

		# check mapping quality
		##INFO=<ID=T_MQ,Number=2,Type=Float,Description="In TUMOR: average mapping quality of consensus indel-supporting reads/reference-supporting reads">
		##INFO=<ID=N_MQ,Number=2,Type=Float,Description="In NORMAL: average mapping quality of consensus indel-supporting reads/reference-supporting reads">
		my ($mq_indel_tum, $mq_ref_tum) = split(",", $x->{INFO}{T_MQ});
		my ($mq_indel_rem, $mq_ref_rem) = split(",", $x->{INFO}{N_MQ});		
		if ($mq_indel_tum < 40)
		{
			INFO("REJECT: POOR MAPPING: ",$x->{CHROM},":",$x->{POS},"\t","T_MQ=$mq_indel_tum,$mq_ref_tum");
			next;
		}
		
				
#		print "reads_indel_fwd: $reads_indel_fwd\n";
#		print "reads_indel_rev: $reads_indel_rev\n";
#		print "reads_ref_fwd: $reads_ref_fwd\n";
#		print "reads_ref_rev: $reads_ref_rev\n";
	}
	else
	{
		croak "ERROR: Invalid variant type: $var_type\n";
	}

	my (@repeats, @dups, @blacklist, @retro, @rem_samples);
	my ($chr, $pos) = ($x->{CHROM}, $x->{POS});

	# ----- annotate overlapping repeat regions
	{
		my $iter = $rmsk->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $rmsk->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@repeats, "$s[10]:$s[11]:$s[12]");
			}
		}		
	}

	{
		my $iter = $simpleRepeat->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $simpleRepeat->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@repeats, "$s[16]($s[6])");
			}
		}		
	}
	$numrep ++ if (@repeats > 0);

	# ----- annotate segmental duplications
	{
		my $iter = $segdup->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $segdup->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@dups, "$s[4]:".sprintf("%.2f", $s[26]));
			}		
		}		
		$numsegdup ++ if (@dups > 0);
	}	
	
	# ----- annotate overlapping DAC blacklisted regions
	{
		my $iter = $blacklistdb->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $blacklistdb->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@blacklist, "$s[4]");
			}		
		}
		$num_blacklist ++ if (@blacklist > 0);		
	}

	# ----- annotate overlapping g1k accessible regions
	my $accessible = "no";
	{
		my $iter = $g1kAcc->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $g1kAcc->read($iter)) 
			{
				$accessible = "";
				last;
			}		
		}
		$num_not_accessible ++ if ($accessible eq "no");		
	}

	# ----- annotate overlapping retrotransposed (pseudo) genes
	{
		my $iter = $ucscRetro->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $ucscRetro->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@retro, $s[10]);
			}		
		}
		$num_retro ++ if (@retro > 0);
	}

	# ----- annotate variants found in remission samples
	{
		my $iter = $remission->query($chr, $pos-1, $pos);
		if ($iter and $iter->{_})
		{
			while (my $line = $remission->read($iter)) 
			{
				my ($sample, $rchr, $rpos, $ref_allele, $alt_allele, $dp, $ad, $gt) = split("\t", $line);
				if ($pos eq $rpos and $x->{REF} eq $ref_allele and $x->{ALT}->[0] eq $alt_allele and $ad >= 3 and $ad/$dp > 0.05)
				{
					push(@rem_samples, "$sample($ad)");
				}
			}		
		}
		$num_remission ++ if (@rem_samples > 0);
	}
	
	my @rejected_because;
	if ($rejected_variants{"$patient\t$cmp_type\t".$x->{CHROM}."\t".$x->{POS}}) { push(@rejected_because, "manual inspection (".$rejected_variants{"$patient\t$cmp_type\t".$x->{CHROM}."\t".$x->{POS}}.")")}
	if (@repeats > 0) { push(@rejected_because, "repetitive region"); }
	if (@dups > 0) { push(@rejected_because, "segmental duplication"); }
	if (@blacklist > 0) { push(@rejected_because, "blacklisted region"); }
	if (@retro > 0) { push(@rejected_because, "retrotransposon"); }
	if (@rem_samples > 1) { push(@rejected_because, "present remissions"); }
	if ($x->{ID} and $x->{ID} ne ".")  { push(@rejected_because, "dbSNP"); }  
	if (defined $x->{INFO}{'dbNSFP_1000Gp1_AF'} and $x->{INFO}{'dbNSFP_1000Gp1_AF'} > 0) { push(@rejected_because, "g1k"); }
	
	$line =~ s/^([^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t)[^\t]+/$1REJECT/ if (@rejected_because > 0);
	
	print VCFOUT "$line" if ($vcf_out);
	
	print "$patient\t";		
	print "$cmp_type\t";
	print "$var_type\t";
	print @rejected_because > 0 ? "REJECT\t" : "$status\t";
	print @rejected_because > 0 ? join(";", @rejected_because) : "", "\t";
	print $x->{CHROM},"\t";
	print $x->{POS},"\t";
	print $x->{ID},"\t";
	print $x->{REF},"\t";
	print $x->{ALT}->[0],"\t";
	my ($gene, $add_genes, $impact, $effect, $affected_exon, $aa_change) = get_impact($x->{INFO}{EFF});
	print "$gene\t";
	print "$add_genes\t";
	print "$impact\t";
	print "$effect\t";
	print "$affected_exon\t";
#	print join(",", @{$x->{FILTER}}),"\t";
#	print exists $x->{gtypes}{$cmp_sample}{SS} ? $variant_stati{$x->{gtypes}{$cmp_sample}{SS}} : "n/a", "\t";
	print "$dp_rem\t";
	print "$ad_rem_ref\t";
	print "$ad_rem_alt\t";
	print "$freq_rem\t";
	print "$dp_tum\t";
	print "$ad_tum_ref\t";
	print "$ad_tum_alt\t";
	print "$freq_tum\t";
	print "$aa_change\t";
	print "EFF=",$x->{INFO}{EFF},"\t";
	print defined $x->{INFO}{'dbNSFP_Polyphen2_HVAR_pred'} ? $x->{INFO}{'dbNSFP_Polyphen2_HVAR_pred'} : "", "\t"; # Polyphen2 prediction based on HumVar, 'D' ('porobably damaging'), 'P' ('possibly damaging') and 'B' ('benign'). Multiple entries separated by ';' 
	print defined $x->{INFO}{'dbNSFP_SIFT_score'} ? $x->{INFO}{'dbNSFP_SIFT_score'} : "", "\t"; # SIFT score, If a score is smaller than 0.05 the corresponding NS is predicted as 'D(amaging)'; otherwise it is predicted as 'T(olerated)'
	print defined $x->{INFO}{'dbNSFP_GERP++_RS'} ? $x->{INFO}{'dbNSFP_GERP++_RS'} : "", "\t"; # GERP++ RS score, the larger the score, the more conserved the site 
	print defined $x->{INFO}{'dbNSFP_29way_logOdds'} ? $x->{INFO}{'dbNSFP_29way_logOdds'} : "", "\t"; # SiPhy score based on 29 mammals genomes. The larger the score, the more conserved the site.
	my $domains = $x->{INFO}{'dbNSFP_Interpro_domain'}; # domain or conserved site on which the variant locates. Domain annotations come from Interpro database. The number in the brackets following a specific domain is the count of times Interpro assigns the variant position to that domain, typically coming from different predicting databases
	if ($domains)
	{
		$domains =~ s/\),/\)\|/g;
		$domains =~ s/\|$//;
		print "$domains\t";
	}
	else
	{
		print "\t";
	}
	print defined $x->{INFO}{'dbNSFP_1000Gp1_AF'} ? $x->{INFO}{'dbNSFP_1000Gp1_AF'} : "";  # Alternative allele frequency in the whole 1000Gp1 data
	print "\t", join(',', @repeats), "\t", join(',', @dups), "\t", join(',', @blacklist), "\t$accessible\t", join(",", @retro), "\t", join(",", @rem_samples);
	print "\n";
		
#	print "\n"; print Dumper($x); exit;
}
$vcf->close();
close(VCFOUT) if ($vcf_out);
	
if ($debug)
{
	INFO("  Total number of variants: $tot_var");
	INFO("  Variants by quality:");
	foreach my $k (keys(%qual_num))
	{
		INFO("    $k: ", $qual_num{$k});
	}
	INFO("  Rejected by MuTect: $filtered_qual");
	INFO("  Excluded due to equal genotype: $filtered_gt");
	INFO("  Excluded due to missing alternative allele: $filtered_alt");
	INFO("  Excluded germline variants: $filtered_germ");
	INFO("  $numrep variants annotated with overlapping repeat.");
	INFO("  $num_blacklist variants annotated with overlapping blacklisted region.");
	INFO("  $numsegdup variants annotated with overlapping segmental duplication.");
	INFO("  $num_not_accessible variants fall into G1K non-accessible region.");
	INFO("  $num_retro variants annotated with overlapping retrotransposed (pseudo)gene.");
	INFO("  $num_remission variants present in remission sample(s).");
}

# ------------------------------------------

sub get_impact
{
	my $effs = shift or die "ERROR: effect not specified";

	# determine all genes impacted by variants
	my (%genes_by_impact, %all_genes, $combined_effect, $combined_impact, %affected_exons, %aa_changes);
	foreach my $eff (split(",", $effs))
	{
		my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
			or die "ERROR: could not parse SNP effect: $effs\n";

		my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
			$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
				or die "ERROR: could not parse SNP effect: $eff\n";
		 
		$aa_changes{$aa_change} = 1 if ($aa_change);

		if ($exon and $transcript and $gene_name)
		{
			$transcript =~ s/\.\d+$//; # remove version number from accession
			$transcript =~ s/\.\d+$//; 
			$affected_exons{$gene_name}{$exon}{$transcript} = 1;
			if ($canonical{$transcript})
			{
				$affected_exons{$gene_name}{'canonical'}{$exon}{$transcript} = 1;
			}
		}
			
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

	my @aff_exons;
	foreach my $g (keys(%affected_exons))
	{
		if (exists $affected_exons{$g}{'canonical'}) # known canonical transcript for this gene?
		{
			foreach my $e (keys(%{$affected_exons{$g}{'canonical'}}))
			{
				my @transcripts;
				foreach my $t (keys(%{$affected_exons{$g}{'canonical'}{$e}}))
				{
					push(@transcripts, "$g:$t");
				}
				push(@aff_exons, "$e (".join(";", @transcripts).")");
			}
		}
		else
		{
			foreach my $e (keys(%{$affected_exons{$g}}))
			{
				next if ($e eq 'canonical');

				my @transcripts;
				foreach my $t (keys(%{$affected_exons{$g}{$e}}))
				{
					push(@transcripts, "$g:$t");
				}
				push(@aff_exons, "$e (".join(";", @transcripts).")");
			}
			
		}
	}

	return ($gene, $add_genes, $combined_impact, $combined_effect, 
			@aff_exons > 0 ? join(",", @aff_exons) : "", join(";", keys(%aa_changes)));
}
