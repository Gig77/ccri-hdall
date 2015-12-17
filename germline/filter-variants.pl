use warnings FATAL => qw( all );
use strict;

use lib "/mnt/projects/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Tabix;
use Vcf;
use Data::Dumper;
use Getopt::Long;
use Carp;

my ($vcf_out, $header, $rejected_variants_file, $patient);
my ($rmsk_file, $simplerepeat_file, $blacklist_file, $segdup_file, $g1k_accessible_file, $cosmic_mutation_file, $ucsc_retro_file, $remission_variants_file, $evs_file, $clinvar_file);
GetOptions
(
	"patient=s" => \$patient,  # patient ID
	"vcf-out=s" => \$vcf_out,  # filtered VCF output file
	"header" => \$header,  # if set, write header line to output
	"rmsk-file=s" => \$rmsk_file, # TABIX indexed UCSC table rmsk
	"simpleRepeat-file=s" => \$simplerepeat_file, # TABIX indexed UCSC table rmsk
	"blacklist-file=s" => \$blacklist_file, # TABIX indexed UCSC table wgEncodeDacMapabilityConsensusExcludable
	"segdup-file=s" => \$segdup_file, # TABIX indexed UCSC table genomicSuperDups
	"g1k-accessible=s" => \$g1k_accessible_file, # TABIX indexed UCSC table tgpPhase1AccessibilityPilotCriteria
	"ucscRetro=s" => \$ucsc_retro_file, # TABIX indexed UCSC table ucscRetroAli5
	"remission-variants-file=s" => \$remission_variants_file, # TABIX indexed file with variants found in remission samples (GATK)
	"cosmic-mutation-file=s" => \$cosmic_mutation_file,
	"evs-file=s" => \$evs_file, # TABIX indexed file with wariants from Exome Variant Server (http://evs.gs.washington.edu/EVS/)
	"clinvar-file=s" => \$clinvar_file # TABIX indexed VCF file with ClinVar variants (ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf)
);

# TABLE: filtered-variants
if ($header)
{
	print "patient\t";		
	print "var_type\t";
	print "type\t";
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
	print "non_silent\t";
	print "deleterious\t";
	print "exons\t";
	print "dp_rem_tot\t";
	print "dp_rem_ref\t";
	print "dp_rem_var\t";
	print "freq_rem\t";
	print "both_strands_rem\t";
	print "pval_rem\t";
	print "dp_leu_tot\t";
	print "dp_leu_ref\t";
	print "dp_leu_var\t";
	print "freq_leu\t";
	print "both_strands_leu\t";
	print "pval_leu\t";
	print "aa_change\t";
	print "effect_long\t";
	print "Polyphen2\t";
	print "SIFT\t";
	print "GERP++\t";
	print "SiPhy\t";
	print "InterPro\t";
	print "AF_1000G\t";
	print "cosmic_hits_nt\t";
	print "cosmic_hits_aa\t";
	print "cosmic_hits_leu_nt\t";
	print "cosmic_hits_leu_aa\t";
	print "repeat\t";
	print "tsegdup\t";
	print "blacklist\t";
	print "g1k-accessible\t";	
	print "evs_variant\t";	
	print "clinvar_variant\n";	
	exit;
}

my $debug = 1;

my $vcf_file = $ARGV[0] or croak "ERROR: VCF file not specified\n";

croak "ERROR: --patient not specified" if (!$patient);
croak "ERROR: --rmsk-file not specified" if (!$rmsk_file);
croak "ERROR: --simpleRepeat-file not specified" if (!$simplerepeat_file);
croak "ERROR: --blacklist-file not specified" if (!$blacklist_file);
croak "ERROR: --segdup-file not specified" if (!$segdup_file);
croak "ERROR: --g1k-accessible not specified" if (!$g1k_accessible_file);
croak "ERROR: --cosmic-mutation-file not specified" if (!$cosmic_mutation_file);
croak "ERROR: --ucscRetro not specified" if (!$ucsc_retro_file);
croak "ERROR: --remission-variants-file not specified" if (!$remission_variants_file);
croak "ERROR: --evs-file not specified" if (!$evs_file);
croak "ERROR: --clinvar-file not specified" if (!$clinvar_file);


my $tum_sample = "TUMOR"; 
my $rem_sample = "NORMAL"; 

# read kgXref, knownCanonical to determine UCSC canonical transcripts affected by variant
my %kgID2refSeq;
open(G,"/mnt/projects/generic/data/hg19/hg19.kgXref.txt") or die "could not open file /mnt/projects/generic/data/hg19/hg19.kgXref.txt";
while(<G>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refSeq, $protAcc, $description, $rfamAcc, $tRnaName) = split(/\t/);

	$kgID2refSeq{$kgID} = $refSeq if ($refSeq);
}
close(G);
INFO(scalar(keys(%kgID2refSeq))." gene descriptions read from file /mnt/projects/generic/data/hg19/hg19.kgXref.txt");

my %canonical;
open(G,"/mnt/projects/generic/data/hg19/hg19.knownCanonical.txt") or die "could not open file /mnt/projects/generic/data/hg19/hg19.knownCanonical.txt";
<G>; # skip header
while(<G>)
{
	chomp;
	my ($chrom, $chromStart, $chromEnd, $clusterId, $transcript, $protein) = split(/\t/);
	
	$canonical{$kgID2refSeq{$transcript}} = 1 if ($kgID2refSeq{$transcript});
}
close(G);
INFO(scalar(keys(%canonical))." canonical genes read from file /mnt/projects/generic/data/hg19/hg19.knownCanonical.txt");

# read cosmic mutations
my (%cosmic, %cosmic_leuk);
my $entries_read = 0;
open(C, "$cosmic_mutation_file") or croak "ERROR: Could not open file $cosmic_mutation_file";
<C>; # skip header
while(<C>)
{
	chomp;
	my ($gene_name, $accession_number, $gene_cds_length, $hgnc_id, $sample_name, $id_sample, $id_tumour, $primary_site, $site_subtype, $primary_histology,
		$histology_subtype, $genome_wide_screen, $mutation_id, $mutation_cds, $mutation_aa, $mutation_description, $mutation_zygosity, $mutation_ncbi36_genome_position,
		$mutation_ncbi36_strand, $mutation_GRCh37_genome_position, $mutation_GRCh37_strand, $mutation_somatic_status, $pubmed_pmid, $sample_source, 
		$tumour_origin, $age, $comments) = split /\t/;
	
	next if ($mutation_somatic_status ne "Confirmed somatic variant");
	$gene_name =~ s/_ENST.*//;
	
	$cosmic{$mutation_GRCh37_genome_position} = defined $cosmic{$mutation_GRCh37_genome_position} ? $cosmic{$mutation_GRCh37_genome_position} + 1 : 1;
	$cosmic_leuk{$mutation_GRCh37_genome_position} = defined $cosmic_leuk{$mutation_GRCh37_genome_position} ? $cosmic_leuk{$mutation_GRCh37_genome_position} + 1 : 1
		if ($histology_subtype =~ /leukaemia/);
	
	if ($mutation_aa =~ /p\.(.)(\d+)(.+)/)
	{
		my ($prev_aa, $aa_pos, $after_aa) = ($1, $2, $3);
		$cosmic{"$gene_name:$prev_aa:$aa_pos"} = defined $cosmic{"$gene_name:$prev_aa:$aa_pos"} ? $cosmic{"$gene_name:$prev_aa:$aa_pos"} + 1 : 1;
		$cosmic_leuk{"$gene_name:$prev_aa:$aa_pos"} = defined $cosmic_leuk{"$gene_name:$prev_aa:$aa_pos"} ? $cosmic_leuk{"$gene_name:$prev_aa:$aa_pos"} + 1 : 1
			if ($histology_subtype =~ /leukaemia/);
	}
	$entries_read ++;
}
close(C);
INFO("$entries_read mutations read from file $cosmic_mutation_file");


my $rmsk = Tabix->new(-data => $rmsk_file);
my $simpleRepeat = Tabix->new(-data => $simplerepeat_file);
my $blacklistdb = Tabix->new(-data => $blacklist_file);
my $segdup = Tabix->new(-data => $segdup_file);
my $g1kAcc = Tabix->new(-data => $g1k_accessible_file);
my $ucscRetro = Tabix->new(-data => $ucsc_retro_file);
my $remission = Tabix->new(-data => $remission_variants_file);
my $evs = Tabix->new(-data => $evs_file);
my $clinvar = Tabix->new(-data => $clinvar_file);

$| = 1; # turn on autoflush

INFO("Processing file $vcf_file...");

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
die "ERROR: Sample name $tum_sample not found!\n" if ($tum_sample ne $samples[0] and $tum_sample ne $samples[1]);
die "ERROR: Sample name $rem_sample not found!\n" if ($rem_sample ne $samples[0] and $rem_sample ne $samples[1]);

my ($tot_var, $filtered_alt, $filtered_germ) = (0, 0, 0, 0, 0);
my ($numrep, $num_blacklist, $numsegdup, $num_not_accessible, $num_dbsnp, $num_retro, $num_remission, $num_evs, $num_clinvar) = (0, 0, 0, 0, 0, 0, 0, 0, 0);
my %qual_num;

my %variant_stati = 
(
	0 => 'reference',
	1 => 'germline',
	2 => 'somatic',
	3 => 'LOH',
	5 => 'unknown'
);

##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Variant Clinical Significance, 0 - unknown, 1 - untested, 2 - non-pathogenic, 3 - probable-non-pathogenic, 4 - probable-pathogenic, 5 - pathogenic, 6 - drug-response, 7 - histocompatibility, 255 - other">
my %clinvar_sig =
(
	0 => "unknown", 
	1 => "untested",
	2 => "non-pathogenic",
	3 => "probable-non-pathogenic",
	4 => "probable-pathogenic",
	5 => "pathogenic",
	6 => "drug-response",
	7 => "histocompatibility",
	255 => "other"
);

while (my $line = $vcf->next_line())
{
	my $x = $vcf->next_data_hash($line);

	$tot_var ++;
	$qual_num{$x->{FILTER}->[0]} = $qual_num{$x->{FILTER}->[0]} ? $qual_num{$x->{FILTER}->[0]} + 1 : 1;
	
	my $gt_tum = $x->{gtypes}{$tum_sample}{GT};
	die "ERROR: Could not determine genotype of sample $tum_sample in file $vcf_file\n" if (!defined $gt_tum or $gt_tum eq "");

	if (@{$x->{ALT}} != 1) # more than one alternative allele?
	{
		$filtered_alt ++;
		next;
	}		

	my $ref_allele = $x->{REF};
	my $alt_allele = $x->{ALT}->[0];
	my $var_type;
	if (length($ref_allele) == length($alt_allele))
	{
		$var_type = 'snp';
	}
	elsif (length($ref_allele) < length($alt_allele))
	{
		$var_type = 'ins';
	}
	else
	{
		$var_type = 'del';
	}

	my $status = $x->{FILTER}->[0];		
	
	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">				
	my $dp_tum = $x->{gtypes}{$tum_sample}{DP};
	my $dp_rem = $x->{gtypes}{$rem_sample}{DP};

	##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
	##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
	my $ad_tum_ref = $x->{gtypes}{$tum_sample}{RD};
	my $ad_tum_alt = $x->{gtypes}{$tum_sample}{AD};
	my $ad_rem_ref = $x->{gtypes}{$rem_sample}{RD};
	my $ad_rem_alt = $x->{gtypes}{$rem_sample}{AD};
	
	##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
	my $freq_tum = $x->{gtypes}{$tum_sample}{FREQ};
	$freq_tum =~ s/\%//;
	my $freq_rem = $x->{gtypes}{$rem_sample}{FREQ};
	$freq_rem =~ s/\%//;

	##INFO=<ID=SPV,Number=1,Type=Float,Description="Fisher's Exact Test P-value of tumor versus normal for Somatic/LOH calls">
	##INFO=<ID=GPV,Number=1,Type=Float,Description="Fisher's Exact Test P-value of tumor+normal versus no variant for Germline calls">
	my $pval_tum = $x->{INFO}{'SPV'};
	my $pval_rem = $x->{INFO}{'GPV'};
			
	##FORMAT=<ID=DP4,Number=1,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev">
	my ($ref_fwd_tum, $ref_rev_tum, $var_fwd_tum, $var_rev_tum) = split(",", $x->{gtypes}{$tum_sample}{DP4});		
	my ($ref_fwd_rem, $ref_rev_rem, $var_fwd_rem, $var_rev_rem) = split(",", $x->{gtypes}{$rem_sample}{DP4});
	
	my $both_strands_tum = ($var_fwd_tum and $var_rev_tum);
	my $both_strands_rem = ($var_fwd_rem and $var_rev_rem);

	next if ($var_fwd_tum < 2 || $var_rev_tum < 2);
	next if ($dp_rem < 10);
	next if (!$both_strands_rem or !$both_strands_tum);

	my (@repeats, @dups, @blacklist, @retro, @rem_samples, %evss, @clinvars);
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

	# ----- annotate variants found in Exome Variant Server
	{
		my $iter = $evs->query($chr, $pos-1, $pos);
		if ($iter and $iter->{_})
		{
			while (my $line = $evs->read($iter)) 
			{
				my ($echr, $epos, $rsID, $dbSNPVersion, $alleles, $europeanAmericanAlleleCount, $africanAmericanAlleleCount, $allAlleleCount, $MAFinPercent_EA_AA_All, $europeanAmericanGenotypeCount, 
					$africanAmericanGenotypeCount, $allGenotypeCount, $avgSampleReadDepth, $genes, $geneAccession, $functionGVS, $hgvsProteinVariant, $hgvsCdnaVariant, $codingDnaSize, 
					$conservationScorePhastCons, $conservationScoreGERP, $granthamScore, $polyphen2_score, $refBaseNCBI37, $chimpAllele, $clinicalInfo, $filterStatus, $onIlluminaHumanExomeChip,
					$gwasPubMedInfo, $EA_EstimatedAge_kyrs, $AA_EstimatedAge_kyrs) = split(/\s/, $line);
					
				next if ($echr ne $chr or $epos ne $pos);
				foreach my $allele (split(";", $alleles))
				{
					my ($ref, $alt) = $allele =~ /(.+)\>(.+)/;
					if ($ref eq $x->{REF} and $alt eq $x->{ALT}->[0])
					{
						my ($alt_count, $ref_count) = $allAlleleCount =~ /(\d+).+?(\d+)/;
						my $alt_percent = sprintf("%.3f", $alt_count/($alt_count+$ref_count)*100);
						$evss{$alt_percent} = 1;					
					}
				}
			}		
		}
		$num_evs ++ if (keys(%evss) > 0);
	}

	# ----- ClinVar variants
	{
		my $schr = $chr =~ s/^chr//;
		my $iter = $clinvar->query($schr, $pos-1, $pos);
		if ($iter and $iter->{_})
		{
			while (my $line = $clinvar->read($iter)) 
			{
				my ($echr, $epos, $rsID, $eref, $ealt) = split(/\s/, $line);
				if ($schr eq $echr and $pos eq $epos and $x->{REF} eq $eref and $x->{ALT}->[0] eq $ealt)
				{
					my ($acc) = $line =~ /CLNACC=([^;\n]+)/;
					my ($sig) = $line =~ /CLNSIG=(\d+)/;
					$sig = $clinvar_sig{$sig};
					push(@clinvars, "$acc(".$sig.")");
				}
			}		
		}
		$num_clinvar ++ if (@clinvars > 0);
	}

	my $loc = $x->{CHROM}.":".$x->{POS}."-".$x->{POS};
	$loc =~ s/^chr//;
	
	my $reject = 0;
	my @reject_because;

	#if ($variant_stati{$x->{INFO}{'SS'}} eq "SOMATIC") { $reject = 1; push(@reject_because, "somatic"); }
	if (@repeats > 0) { $reject = 1; push(@reject_because, "repetitive region"); };
	if (@dups > 0) { $reject = 1; push(@reject_because, "segmental duplication"); };
	if (@blacklist > 0) { $reject = 1; push(@reject_because, "blacklisted region"); };
	#if (@retro > 0) { $reject = 1; push(@reject_because, "retrotransposon"); }
	#if (@rem_samples > 1) { $reject = 1; push(@reject_because, "present remissions"); }
	if ($x->{ID} and $x->{ID} ne ".")  { $num_dbsnp ++; };
	
	next if ($reject);  

	my $polyphen = $x->{INFO}{'dbNSFP_Polyphen2_HVAR_pred'};
	my $sift = $x->{INFO}{'dbNSFP_SIFT_score'};
	my $siphy = $x->{INFO}{'dbNSFP_29way_logOdds'};
	my ($gene, $add_genes, $impact, $effect, $affected_exon, $aa_change) = get_impact($x->{INFO}{EFF});
	
	my $is_deleterious = "n/d";
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and $polyphen and $polyphen =~ /D/ and defined $sift and $sift < 0.05); # polyphen & sift
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and $polyphen and $polyphen =~ /D/ and defined $siphy and $siphy >= 12); # polyphen & siphy
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and defined $sift and $sift < 0.05 and defined $siphy and $siphy >= 12); # sift and siphy
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and defined $siphy and $siphy > 20); # siphy only, highly conserved (keeps GNAQ)
	$is_deleterious = "yes" if ($effect eq "NON_SYNONYMOUS_CODING" and defined $siphy and $siphy > 20); # siphy only, highly conserved (keeps GNAQ)
	$is_deleterious = "yes" if ($effect eq "FRAME_SHIFT" or $effect eq "SPLICE_SITE_ACCEPTOR" or $effect eq "SPLICE_SITE_DONOR" or $effect eq "STOP_GAINED");
	$is_deleterious = "no" if ($is_deleterious ne "yes" and defined $polyphen and defined $sift);
	$is_deleterious = "no" if ($effect eq "DOWNSTREAM" or $effect eq "UPSTREAM" or $effect eq "INTRON" or $effect eq "INTERGENIC" or $effect eq "SYNONYMOUS_CODING" or $effect eq "SYNONYMOUS_STOP" or $effect eq "SYNONYMOUS_START" or $effect eq "UTR_3_PRIME" or $effect eq "UTR_5_PRIME" or $effect eq "UTR_5_DELETED" or $effect eq "UTR_3_DELETED" or $effect eq "START_GAINED");

	my $non_silent = 0;
	$non_silent = 1 if ($effect eq "STOP_GAINED" or $effect eq "STOP_LOST" or $effect eq "SPLICE_SITE_DONOR" or $effect eq "SPLICE_SITE_ACCEPTOR" or $effect eq "FRAME_SHIFT" or $effect eq "CODON_CHANGE_PLUS_CODON_INSERTION" or $effect eq "CODON_DELETION" or $effect eq "NON_SYNONYMOUS_CODING" or $effect eq "CODON_INSERTION" or $effect eq "CODON_CHANGE_PLUS_CODON_DELETION" or $effect eq "NON_SYNONYMOUS_START" or $effect eq "START_LOST");
	#next if (!$non_silent);
	
	$line =~ s/^([^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t)[^\t]+/$1REJECT/ if ($reject);
	
	print VCFOUT "$line" if ($vcf_out);
	
	print "$patient\t";		
	print "$var_type\t";
	print $variant_stati{$x->{INFO}{'SS'}}, "\t";
	print $reject ? "REJECT\t" : "$status\t";
	print join(";", @reject_because), "\t";
	print $x->{CHROM},"\t";
	print $x->{POS},"\t";
	print $x->{ID},"\t";
	print $x->{REF},"\t";
	print $x->{ALT}->[0],"\t";

	print "$gene\t";
	print "$add_genes\t";
	print "$impact\t";
	print "$effect\t";
	print "$non_silent\t";
	print "$is_deleterious\t";
	print "$affected_exon\t";

	print "$dp_rem\t";
	print "$ad_rem_ref\t";
	print "$ad_rem_alt\t";
	print "$freq_rem\t";
	print $both_strands_rem ? "yes" : "no", "\t";
	print "$pval_rem\t";

	print "$dp_tum\t";
	print "$ad_tum_ref\t";
	print "$ad_tum_alt\t";
	print "$freq_tum\t";
	print $both_strands_tum ? "yes" : "no", "\t";
	print "$pval_tum\t";
	
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
	print defined $x->{INFO}{'dbNSFP_1000Gp1_AF'} ? $x->{INFO}{'dbNSFP_1000Gp1_AF'} : "", "\t";  # Alternative allele frequency in the whole 1000Gp1 data
	print $cosmic{$loc} ? $cosmic{$loc} : "0", "\t";
	print aa_hits([$gene, split(",", $add_genes)], "EFF=".$x->{INFO}{EFF}), "\t";	
	print $cosmic_leuk{$loc} ? $cosmic_leuk{$loc} : "0", "\t";
	print aa_hits([$gene, split(",", $add_genes)], "EFF=".$x->{INFO}{EFF}, 1);	
	
	print "\t", join(',', @repeats), "\t", join(',', @dups), "\t", join(',', @blacklist), "\t$accessible", "\t", join(";", keys(%evss));
	print "\t", join(";", @clinvars);
	print "\n";		
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
	INFO("  Excluded germline variants: $filtered_germ");
	INFO("  Excluded due to missing alternative allele: $filtered_alt");
	INFO("  $numrep variants annotated with overlapping repeat.");
	INFO("  $num_blacklist variants annotated with overlapping blacklisted region.");
	INFO("  $numsegdup variants annotated with overlapping segmental duplication.");
	INFO("  $num_not_accessible variants fall into G1K non-accessible region.");
	INFO("  $num_dbsnp variants common non-pathogenic dbSNP variants.");
	INFO("  $num_retro variants annotated with overlapping retrotransposed (pseudo)gene.");
	INFO("  $num_remission variants present in remission sample(s).");
	INFO("  $num_evs variants present in Exome Variant Server (EVS).");
	INFO("  $num_clinvar variants present in ClinVar.");
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

sub aa_hits
{
	my $genes = shift;
	my $snpeff = shift;
	my $leuk = shift;
	
	return "non-coding" if (@$genes == 0);
	
	my $aa_change_found = 0; 
	foreach my $gene (@$genes) # check each gene
	{
		foreach my $eff (split(",", $snpeff)) # check all isoforms for cosmic match
		{
			my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
				or croak "ERROR: could not parse SNP effect: $snpeff";
	
			my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
				$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
					or croak "ERROR: could not parse SNP effect: $eff";
					 
			if ($aa_change =~ /(.)(\d+)(.+)/)
			{
				$aa_change_found = 1;			
				my ($prev_aa, $aa_pos, $after_aa) = ($1, $2, $3);
				return $cosmic{"$gene:$prev_aa:$aa_pos"} if (!$leuk and defined $cosmic{"$gene:$prev_aa:$aa_pos"});
				return $cosmic_leuk{"$gene:$prev_aa:$aa_pos"} if ($leuk and defined $cosmic_leuk{"$gene:$prev_aa:$aa_pos"});
			}
		}				
	}
	
	return "non-coding" if (!$aa_change_found);
	return "0";
}