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

my ($music_roi, $keep_variants_file, $entrez_mapping, $sample, $header);
GetOptions
(
	"music-roi=s" => \$music_roi,  # MuSiC region of interest (ROI) file (must be tabix accessible, i.e. compressed and indexed)
#	"keep-variants-file=s" => \$keep_variants_file,  # tab-separated file with variants to keep (chr, start)
	"mapping-entrez=s" => \$entrez_mapping,  # file with mappings from gene symbol to entrez ids
	"sample=s" => \$sample,  # e.g. 314_rem_dia
	"header" => \$header  # output header yes/no
);

if ($header)
{
	print_header();
	exit;
}

die "ERROR: --music-roi not specified (.gz file)\n" if (!$music_roi);
die "ERROR: --mapping-entrez not specified\n" if (!$entrez_mapping);
die "ERROR: --sample not specified\n" if (!$sample);

my ($patient, $sample_normal, $sample_tumor) = split("_", $sample) or die "ERROR: invalid sample\n";
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
		$previous_symbols, $previous_names, $aliases, $name_aliases, $entrez_gene_id_ncbi) = split("\t");

	$entrez_gene_id = $entrez_gene_id_ncbi if (!$entrez_gene_id);
	next if (!$entrez_gene_id);
	
	$sym2entrez{$approved_symbol} = $entrez_gene_id;
	map {$sym2entrez{$_} = $entrez_gene_id } split(", ", $previous_symbols);
}
close(GENES);

my $vcf = Vcf->new(file => "-");
$vcf->parse_header();
my (@samples) = $vcf->get_samples();

while (my $x = $vcf->next_data_hash())
{
	my $chr = $x->{CHROM};
    $chr =~ s/^chr//;
	my $pos = $x->{POS};
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

	INFO("$sample_tumor: Variant $chr:$pos impacting multiple genes: ".join(",", keys(%$snpeff_genes)))
		if (keys(%$snpeff_genes) > 1);

	my $rois = get_rois($chr, $pos, $pos+1);

	INFO("$sample_tumor: Variant $chr:$pos mapping to multiple ROIs: ".join(",", keys(%$rois)))
		if (keys(%$rois) > 1);
		
#	my %written;
	foreach my $gene (keys(%$snpeff_genes))
	{
#		next if (exists $written{$gene});
		
		my $effect = get_variant_classification($snpeff_genes->{$gene}, $var_type);
		
		if ($effect eq 'Intron')
		{
			INFO("$sample_tumor: Variant $chr:$pos NOT written: mapping to intron [SnpEff=$gene(".$snpeff_genes->{$gene}, ")]");
			next;
		}
		
		if (!exists $rois->{$gene})
		{
			#my $roi = (keys(%$rois))[0];
			INFO("$sample_tumor: Variant $chr:$pos NOT written: not mapping to ROI [SnpEff=$gene(".$snpeff_genes->{$gene}, ")]");
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
		print "\n";
		
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
	print "Matched_Norm_Sample_UUID\n";	
}

sub get_impacted_genes
{
	my $effs = shift or die "ERROR: effect not specified";

	my %genes;
	foreach my $eff (split(",", $effs))
	{
		my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
			or croak "ERROR: could not parse SNP effect: $effs";

		my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
			$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
				or croak "ERROR: could not parse SNP effect: $eff"; 

		$genes{$gene_name} = $effect
			if ($gene_name);
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