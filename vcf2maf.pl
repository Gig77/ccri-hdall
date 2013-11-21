use warnings FATAL => qw( all );
use strict;
use Carp;

use lib "$ENV{HOME}/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Vcf;
use Tabix;
use Getopt::Long;
use Data::Dumper;

#Log::Log4perl->get_logger()->level($ERROR);

INFO("START: ", join(" ", @ARGV));

my ($music_roi, $keep_variants_file, $entrez_mapping, $sample, $header, $min_af, $deleterious);
GetOptions
(
	"sample=s" => \$sample,  # e.g. 314_rem_dia
	"music-roi=s" => \$music_roi,  # MuSiC region of interest (ROI) file (must be tabix accessible, i.e. compressed and indexed)
#	"keep-variants-file=s" => \$keep_variants_file,  # tab-separated file with variants to keep (chr, start)
	"mapping-entrez=s" => \$entrez_mapping,  # file with mappings from gene symbol to entrez ids
	"header" => \$header,  # output header yes/no
	"deleterious" => \$deleterious,  # output only variants predicted to be deleterious by PolyPhen or SIFT
	"min-af=s" => \$min_af  # minimum allelic frequency
);

if ($header)
{
	print_header();
	exit;
}

die "ERROR: --sample not specified\n" if (!$sample);
die "ERROR: --music-roi not specified (.gz file)\n" if (!$music_roi);
die "ERROR: --mapping-entrez not specified\n" if (!$entrez_mapping);
die "ERROR: ROI file does not exist: $music_roi\n" if (!-e $music_roi);

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

my ($patient, $sample_normal, $sample_tumor) = split("_", $sample) or die "ERROR: invalid sample\n";
my $sample_normal_vcf = $sample_normal; 
my $sample_tumor_vcf = $sample_tumor; 
$sample_normal = $patient."_$sample_normal";
$sample_tumor = $patient."_$sample_tumor";

my $roi = Tabix->new(-data => "$music_roi");


# read biomart id mapping
my %sym2entrez;
open(GENES, $entrez_mapping) or die "ERROR: could not read gene list\n";
<GENES>;
while(<GENES>)
{
	chomp;
	my ($approved_symbol, $entrez_gene_id, $accession_numbers, $approved_name, 
		$previous_symbols, $previous_names, $aliases, $name_aliases, $entrez_gene_id_ncbi,
		$ensembl_id, $uniprot_id, $refseq_id, $ucsc_id) = split("\t");

	$entrez_gene_id = $entrez_gene_id_ncbi if (!$entrez_gene_id);
	next if (!$entrez_gene_id);
	
	$sym2entrez{$approved_symbol} = $entrez_gene_id;
	map {$sym2entrez{$_} = $entrez_gene_id } split(", ", $previous_symbols);
}
close(GENES);

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

my (%canonical, %sym2size);
open(G,"$ENV{HOME}/generic/data/hg19/hg19.knownCanonical.txt") or die "could not open file $ENV{HOME}/generic/data/hg19/hg19.knownCanonical.txt";
while(<G>)
{
	chomp;
	my ($chrom, $chromStart, $chromEnd, $clusterId, $transcript, $protein) = split(/\t/);
	
	my $geneSymbol = $id2sym{$transcript};
	my $size = $chromEnd-$chromStart;
	
	# if multiple canonical transcripts for this gene symbol, use larger one
	next if ($canonical{$geneSymbol} and $sym2size{$geneSymbol} > $size); 
	
	$canonical{$geneSymbol} = $transcript;
	$sym2size{$geneSymbol} = $size;
}
close(G);

# next kgXref
my %refseq2ucsc;
open(G,"$ENV{HOME}/generic/data/hg19/hg19.kgXref.txt") or die "could not open file $ENV{HOME}/generic/data/hg19/hg19.kgXref.txt";
while(<G>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refSeq, $protAcc, $description, $rfamAcc, $tRnaName) = split(/\t/);
	next if ($refseq2ucsc{$refSeq} and $canonical{$geneSymbol} and $canonical{$geneSymbol} ne $kgID); # only map canonical transcripts
	$refseq2ucsc{$refSeq} = $kgID;
}
close(G);

my $vcf = Vcf->new(file => "-");
$vcf->parse_header();

my $mutect = $vcf->get_header_line(key => 'GATKCommandLine', ID => 'MuTect')->[0]->{'CommandLineOptions'};
$mutect = $vcf->get_header_line(key => 'MuTect')->[0]->[0]->{'value'} if (!$mutect);
if ($mutect)
{
	($sample_normal_vcf) = $mutect =~ /normal_sample_name=(\S+)/;
	($sample_tumor_vcf) = $mutect =~ /tumor_sample_name=(\S+)/;
}

my (@samples) = $vcf->get_samples();

if ($samples[0] =~ /Diagnosis/ or $samples[1] =~ /Diagnosis/) { $sample_tumor_vcf =~ s/dia/Diagnosis/; }
if ($samples[0] =~ /Relapse/ or $samples[1] =~ /Relapse/) { $sample_tumor_vcf =~ s/rel/Relapse/; }
if ($samples[0] =~ /Remission/ or $samples[1] =~ /Remission/) { $sample_normal_vcf =~ s/rem/Remission/; }

if ($samples[0] =~ /1020540_Diagosnosis/ or $samples[1] =~ /1020540_Diagosnosis/) { $sample_tumor_vcf = "1020540_Diagosnosis"; } # typo
if ($samples[0] =~ /G_Diagnosis_/ or $samples[1] =~ /G_Diagnosis_/) { $sample_tumor_vcf = "G_Diagnosis_"; } # typo
if ($samples[0] =~ /K_Diagnosis_/ or $samples[1] =~ /K_Diagnosis_/) { $sample_tumor_vcf = "K_Diagnosis_"; } # typo
if ($samples[0] =~ /715_Relapse_2/ or $samples[1] =~ /715_Relapse_2/) { $sample_tumor_vcf = "715_Relapse_2"; }

if ($sample_normal_vcf ne $samples[0] and $sample_normal_vcf ne $samples[1])
{
	$sample_normal_vcf = $patient2sample{$patient."_$sample_normal_vcf"} ? $patient2sample{$patient."_$sample_normal_vcf"} : $patient."_$sample_normal_vcf"; 
}
if ($sample_tumor_vcf ne $samples[0] and $sample_tumor_vcf ne $samples[1])
{
	$sample_tumor_vcf = $patient2sample{$patient."_$sample_tumor_vcf"} ? $patient2sample{$patient."_$sample_tumor_vcf"} : $patient."_$sample_tumor_vcf"; 
}

print STDERR "Normal sample name: $sample_normal_vcf\n";
print STDERR "Tumor sample name: $sample_tumor_vcf\n";

# sanity checks
die "ERROR: Sample name $sample_normal_vcf not found!\n" if ($sample_normal_vcf ne $samples[0] and $sample_normal_vcf ne $samples[1]);
die "ERROR: Sample name $sample_tumor_vcf not found!\n" if ($sample_tumor_vcf ne $samples[0] and $sample_tumor_vcf ne $samples[1]);
die "ERROR: Sample names identical: $sample_tumor_vcf!\n" if ($sample_tumor_vcf eq $sample_normal_vcf);

while (my $x = $vcf->next_data_hash())
{	
	my $chr = $x->{CHROM};
    $chr =~ s/^chr//;
	my $pos = $x->{POS};

	next if ($x->{FILTER}->[0] eq "REJECT");
	
	my $ref_allele = $x->{REF};
	my $alt_allele_1 = $ref_allele; # lets assume for now that all variants are heterozygous
	my $alt_allele_2 = $x->{ALT}->[0]; # @{$x->{ALT}} > 1 ? $x->{ALT}->[1] : '';
	my $dbSNP = ($x->{ID} and $x->{ID} ne '.') ? $x->{ID} : "";

	my $snpeff_genes = get_impacted_genes($x->{INFO}{EFF});

	my $var_type;
	if (length($ref_allele) == length($alt_allele_2))
	{
		$var_type = 'snp';
	}
	elsif (length($ref_allele) < length($alt_allele_2))
	{
		$var_type = 'ins';
	}
	else
	{
		$var_type = 'del';
	}
	
	if (keys(%$snpeff_genes) == 0)
	{
		INFO("$sample_tumor: Variant $chr:$pos NOT written: not impacting gene.");
		next;
	}

	my $af;
	if ($var_type eq 'snp')
	{
		$af = $x->{gtypes}{$sample_tumor_vcf}{FA};
	}
	else # indel
	{
		##INFO=<ID=T_DP,Number=1,Type=Integer,Description="In TUMOR: total coverage at the site">
		##INFO=<ID=N_DP,Number=1,Type=Integer,Description="In NORMAL: total coverage at the site">
		my ($dp_tum, $dp_rem) = ($x->{INFO}{T_DP}, $x->{INFO}{N_DP});
		
		##INFO=<ID=T_AC,Number=2,Type=Integer,Description="In TUMOR: # of reads supporting consensus indel/any indel at the site">
		##INFO=<ID=N_AC,Number=2,Type=Integer,Description="In NORMAL: # of reads supporting consensus indel/any indel at the site">
		my ($ad_tum_alt, $ad_tum_any_indel) = split(",", $x->{INFO}{T_AC}); 		
		$af = $ad_tum_alt / $dp_tum;
	}
	
	if ($min_af and $af < $min_af)
	{
		INFO("$sample_tumor: Skipped variant $chr:$pos: allelic frequency $af < $min_af");
		next;
	}
	
	INFO("$sample_tumor: Variant $chr:$pos impacting multiple genes: ".join(",", keys(%$snpeff_genes)))
		if (keys(%$snpeff_genes) > 1);

	print STDERR "$chr:$pos\n";
	my $rois = get_rois($chr, $pos, $pos+1);

	INFO("$sample_tumor: Variant $chr:$pos mapping to multiple ROIs: ".join(",", keys(%$rois)))
		if (keys(%$rois) > 1);
		
#	my %written;
	foreach my $gene (keys(%$snpeff_genes))
	{
#		next if (exists $written{$gene});
		
		my $effect = get_variant_classification($snpeff_genes->{$gene}{'effect'}, $var_type);
		
		if ($deleterious and $effect eq "Missense_Mutation")
		{
			my ($polyphen, $sift, $siphy) = ($x->{INFO}{'dbNSFP_Polyphen2_HVAR_pred'}, $x->{INFO}{'dbNSFP_SIFT_score'}, $x->{INFO}{'dbNSFP_29way_logOdds'});
			
			my $is_deletrious = 0;
			$is_deletrious = 1 if ($polyphen and $polyphen =~ /D/ and defined $sift and $sift < 0.05); # polyphen & sift
			$is_deletrious = 1 if ($polyphen and $polyphen =~ /D/ and defined $siphy and $siphy >= 12); # polyphen & siphy
			$is_deletrious = 1 if (defined $sift and $sift < 0.05 and defined $siphy and $siphy >= 12); # sift and siphy
			$is_deletrious = 1 if (defined $siphy and $siphy > 20); # siphy only, highly conserved (keeps GNAQ)
			
			next if (!$is_deletrious);
		}

		
		if ($effect eq 'Intron')
		{
			INFO("$sample_tumor: Variant $chr:$pos NOT written: mapping to intron [SnpEff=$gene(".$snpeff_genes->{$gene}{'effect'}, ")]");
			next;
		}
		
		if (!exists $rois->{$gene})
		{
			#my $roi = (keys(%$rois))[0];
			INFO("$sample_tumor: Variant $chr:$pos NOT written: not mapping to ROI [SnpEff=$gene(".$snpeff_genes->{$gene}{'effect'}, ")]");
			next;
		}

		my $entrez_id = $sym2entrez{$gene};
			
		print $gene. "\t"; #1 Hugo_Symbol
		print $entrez_id ? $entrez_id : "0","\t"; #2 Entrez_Gene_Id
		print "CeMM\t"; #3 Center
		print "hg19\t"; #4  NCBI_Build
		print "$chr\t"; #5 Chromosome
		print "$pos\t"; #6 Start_Position
		print "$pos\t"; #7 End_Position
		print "+\t"; #8 Strand
		print "$effect\t"; #9 Variant_Classification
		print uc($var_type)."\t"; #10 Variant_Type
		print "".($var_type eq "ins" ? "-" : $ref_allele)."\t"; #11 Reference_Allele
		print "".($var_type eq "ins" ? "-" : $ref_allele)."\t"; #12 Tumor_Seq_Allele1
		print "".($var_type eq "del" ? "-" : $alt_allele_2)."\t"; #13 Tumor_Seq_Allele2
		print "$dbSNP\t"; #14 dbSNP_RS
		print "\t"; #15 dbSNP_Val_Status
		print "$sample_tumor\t"; #16 Tumor_Sample_Barcode
		print "$sample_normal\t"; #17 Matched_Norm_Sample_Barcode
		print "$ref_allele\t"; #18 Match_Norm_Seq_Allele1
		print "$ref_allele\t"; #19 Match_Norm_Seq_Allele2
		print "\t"; #20 Tumor_Validation_Allele1
		print "\t"; #21 Tumor_Validation_Allele2
		print "$ref_allele\t"; #22 Match_Norm_Validation_Allele1
		print "$ref_allele\t"; #23 Match_Norm_Validation_Allele2
		print "Unknown\t"; #24 Verification_Status
		print "Unknown\t"; #25 Validation_Status
		print "Somatic\t"; #26 Mutation_Status
		print "\t"; #27 Sequencing_Phase
		print "Capture\t"; #28 Sequence_Source
		print "\t"; #29 Validation_Method
		print "\t"; #30 Score
		print "\t"; #31 BAM_File
		print "Illumina HiSeq\t"; #32 Sequencer
		print "$sample_tumor\t"; #33 Tumor_Sample_UUID
		print "$sample_normal\t"; #34 Matched_Norm_Sample_UUID
		
		print $snpeff_genes->{$gene}{'transcript'} ? $snpeff_genes->{$gene}{'transcript'} : "", "\t";
		print $snpeff_genes->{$gene}{'aa_change'} ? $snpeff_genes->{$gene}{'aa_change'} : "", "\t";
		print $snpeff_genes->{$gene}{'nt_pos'} ? $snpeff_genes->{$gene}{'nt_pos'} : "", "\n";
		
#		$written{$gene} = 1;
	}

}
$vcf->close();

INFO("END");

## ------------------------------------------

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
	print "Matched_Norm_Sample_UUID\t";
	
	# append additional columns required by 'genome music proximity'
	print "transcript_name\t"; # the transcript name, such as NM_000028
	print "amino_acid_change\t"; #  the amino acid change, such as p.R290H
	print "c_position\n"; # the nucleotide position changed, such as c.869
}

sub get_impacted_genes
{
	my $effs = shift or die "ERROR: effect not specified";

	my (%genes, %tmp);
	foreach my $eff (split(",", $effs))
	{
		my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
			or croak "ERROR: could not parse SNP effect: $effs";

		my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
			$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
				or croak "ERROR: could not parse SNP effect: $eff"; 

		next if (!$gene_name);
		
		$transcript =~ s/\.\d+$//;
		
		# figure out nucleotide coordinate of this variant
		my ($aa_pos) = $aa_change =~ /(\d+)/;
		$codon_change =~ /[ATGC]/; # find first uppercase character
		my $codon_pos = $-[0];
		my $nt_pos = ($aa_pos-1)*3+$codon_pos+1 if (defined $aa_pos and defined $codon_pos);
		
		# remember for this gene
		$tmp{$gene_name}{'effect'} = $effect;
		$tmp{$gene_name}{'transcript'} = $transcript;
		$tmp{$gene_name}{'aa_change'} = "p.$aa_change" if ($aa_change);
		$tmp{$gene_name}{'nt_pos'} = "c.$nt_pos" if ($nt_pos);
				
		# remember also for canonical transcript
		if ($transcript and $refseq2ucsc{$transcript} and $canonical{$gene_name} and $canonical{$gene_name} eq $refseq2ucsc{$transcript})
		{
			$genes{$gene_name}{'effect'} = $effect;
			$genes{$gene_name}{'transcript'} = $transcript;
			$genes{$gene_name}{'aa_change'} = "p.$aa_change" if ($aa_change);
			$genes{$gene_name}{'nt_pos'} = "c.$nt_pos" if ($nt_pos);			
		}
	}

	# return effect of last transcript if no canonical transcript found
	foreach my $g (keys(%tmp))
	{
		next if ($genes{$g}{'effect'} or !$tmp{$g}{'effect'});
		
		$genes{$g}{'effect'} = $tmp{$g}{'effect'};
		$genes{$g}{'transcript'} = $tmp{$g}{'transcript'};
		$genes{$g}{'aa_change'} = $tmp{$g}{'aa_change'};
		$genes{$g}{'nt_pos'} = $tmp{$g}{'nt_pos'};
	}

	return \%genes;
}

# map variant to gene using MuSiC ROI file and tabix
sub get_rois
{
	my ($chr, $start, $end) = @_;
	croak "ERROR: bad input parameters" if (!$chr or !$start or !$end);

	my %rois;

	my $iter = $roi->query($chr, $start, $end);
	return \%rois if (!$iter or !$iter->{_});

	while (my $line = $roi->read($iter)) 
	{
		my $roi = (split("\t", $line))[3];  
		$rois{$roi} = 1; 
	}

	return \%rois;
}

sub get_variant_classification
{
	my ($snpEff_effect, $var_type) = @_;

	croak "ERROR: bad input parameters" if (!$snpEff_effect or !$var_type);
	croak "ERROR: invalid variant type: $var_type" if ($var_type !~ /(snp|ins|del)/);
	
	my %translate =
	(
		RARE_AMINO_ACID 		=> 'Missense_Mutation', # The variant hits a rare amino acid thus is likely to produce protein loss of function
		NON_SYNONYMOUS_START	=> 'Missense_Mutation', # Variant causes start codon to be mutated into another start codon (the new codon produces a different AA). 	Atg/Ctg, M/L (ATG and CTG can be START codons)
		START_LOST 				=> 'Missense_Mutation', # Variant causes start codon to be mutated into a non-start codon. 	aTg/aGg, M/R
		NON_SYNONYMOUS_CODING 	=> 'Missense_Mutation', # Variant causes a codon that produces a different amino acid 	Tgg/Cgg, W/R

		SYNONYMOUS_START 	=> 'Silent', # Variant causes start codon to be mutated into another start codon. 	Ttg/Ctg, L/L (TTG and CTG can be START codons)
		SYNONYMOUS_CODING 	=> 'Silent', # Variant causes a codon that produces the same amino acid 	Ttg/Ctg, L/L
		SYNONYMOUS_STOP 	=> 'Silent', # Variant causes stop codon to be mutated into another stop codon. 	taA/taG, */*

		SPLICE_SITE_ACCEPTOR 	=> 'Splice_Site', # The variant hits a splice acceptor site (defined as two bases before exon start, except for the first exon).
		SPLICE_SITE_DONOR 		=> 'Splice_Site', # The variant hits a Splice donor site (defined as two bases after coding exon end, except for the last exon).

		FRAME_SHIFT => $var_type eq 'del' ? 'Frame_Shift_Del' : 'Frame_Shift_Ins', # Insertion or deletion causes a frame shift 	An indel size is not multple of 3

		CODON_INSERTION 					=> 'In_Frame_Ins', # One or many codons are inserted 	An insert multiple of three in a codon boundary
		CODON_CHANGE_PLUS_CODON_INSERTION  	=> 'In_Frame_Ins', # One codon is changed and one or many codons are inserted 	An insert of size multiple of three, not at codon boundary

		CODON_DELETION 						=> 'In_Frame_Del', #	One or many codons are deleted 	A deletion multiple of three at codon boundary
		CODON_CHANGE_PLUS_CODON_DELETION 	=> 'In_Frame_Del', # One codon is changed and one or more codons are deleted 	A deletion of size multiple of three, not at codon boundary

		UPSTREAM 	=> "5'Flank", # Upstream of a gene (default length: 5K bases)
		DOWNSTREAM 	=> "3'Flank", # Downstream of a gene (default length: 5K bases)

		UTR_5_PRIME 	=> "5'UTR", # Variant hits 5'UTR region
		UTR_5_DELETED 	=> "5'UTR", # The variant deletes an exon which is in the 5'UTR of the transcript
		START_GAINED 	=> "5'UTR", # A variant in 5'UTR region produces a three base sequence that can be a START codon.

		UTR_3_PRIME 	=> "3'UTR", # Variant hits 3'UTR region
		UTR_3_DELETED 	=> "3'UTR", # The variant deletes an exon which is in the 3'UTR of the transcript

		STOP_GAINED 	=> 'Nonsense_Mutation', # Variant causes a STOP codon 	Cag/Tag, Q/*
		STOP_LOST 		=> 'Nonstop_Mutation', # Variant causes stop codon to be mutated into a non-stop codon 	Tga/Cga, */R
		
		INTRON 				=> 'Intron', # Variant hits and intron. Technically, hits no exon in the transcript.
		INTRON_CONSERVED 	=> 'Intron', # The variant is in a highly conserved intronic region

		INTERGENIC 				=> 'IGR', # The variant is in an intergenic region
		INTERGENIC_CONSERVED 	=> 'IGR', # The variant is in a highly conserved intergenic region
		INTRAGENIC 				=> 'IGR', # The variant hits a gene, but no transcripts within the gene

		EXON => 'RNA', # The vairant hits an exon.

#		not mapped:		

		CDS => undef, # The variant hits a CDS.
		GENE => undef, # The variant hits a gene.
		TRANSCRIPT => undef, # The variant hits a transcript.
		EXON_DELETED => undef, # A deletion removes the whole exon.
		CODON_CHANGE => undef, # One or many codons are changed 	An MNP of size multiple of 3

#		? => 'Translation_Start_Site',
#		? => 'RNA',
#		? => 'Targeted_Region',
#		? => 'De_novo_Start_InFrame',
#		? => 'De_novo_Start_OutOfFrame',
	);

	my $maf_effect = $translate{$snpEff_effect};

	if (!$maf_effect)
	{
		ERROR("Could not map snpEff effect to MAF effect: $snpEff_effect");
		$maf_effect = $snpEff_effect;
	}
		
	return $maf_effect;
}