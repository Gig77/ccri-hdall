use warnings FATAL => qw( all );
use strict;

use lib "$ENV{HOME}/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Data::Dumper;
use Carp;

my $debug = 1;
$| = 1; # turn on autoflush


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
INFO(scalar(keys(%id2sym))." id mappgins read from file $ENV{HOME}/hdall/results/id-mappings.tsv");

# read gene description
my %sym2info;
open(G,"$ENV{HOME}/generic/data/hg19/hg19.kgXref.txt") or die "could not open file $ENV{HOME}/generic/data/hg19/hg19.kgXref.txt";
while(<G>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refSeq, $protAcc, $description, $rfamAcc, $tRnaName) = split(/\t/);

	my $sym = $id2sym{$kgID};
	$sym2info{$sym}{'description'} = $description;	
}
close(G);
INFO(scalar(keys(%sym2info))." gene descriptions read from file $ENV{HOME}/generic/data/hg19/hg19.kgXref.txt");

# read cosmic
my $cosmic_file = "$ENV{HOME}/generic/data/cosmic/v65/cancer_gene_census.tsv"; 
my %cosmic;
open(C, $cosmic_file) or die "could not open file $cosmic_file\n";
<C>; # skip header
while(<C>)
{
	chomp;
	my ($symbol, $name, $gene_id, $chr, $chr_band, $cancer_somatic_mut, $cancer_germline_mut, $tumour_types_somatic, $tumour_types_germline, 
	$cancer_syndrome, $tissue_type, $cancer_molecular_genetic, $mutation_type, $translocation_partner, $other_germline_mut, $other_syndrome_disease) = split (/\t/);

	my $tumour_types = $tumour_types_somatic;
	$tumour_types .= ", " if ($tumour_types_somatic and $tumour_types_germline);
	$tumour_types .= $tumour_types_germline;
	
	$cosmic{$symbol} = $tumour_types;
}
close(C);
INFO(scalar(keys(%cosmic))." cancer census genes read from file $cosmic_file");

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
INFO(scalar(keys(%canonical))." canonical genes read from file $ENV{HOME}/generic/data/hg19/hg19.knownCanonical.txt");

my $lines = 0;
open(G,"$ENV{HOME}/generic/data/hg19/hg19.knownGene.txt") or die "could not open file $ENV{HOME}/generic/data/hg19/hg19.knownGene.txt";
while(<G>)
{
	chomp;
	my ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $proteinID, $alignID) = split(/\t/);

	$lines++;
	my $geneSymbol = $id2sym{$name};
	
	next if (exists $canonical{$geneSymbol} and $canonical{$geneSymbol} ne $name and exists $sym2info{$geneSymbol}{'cdslen'}); # prefer canonical transcript (if available)
		
	#$sym2info{$prev2sym{$geneSymbol}}{'exons'} = $exonCount if ($prev2sym{$geneSymbol});
	$sym2info{$geneSymbol}{'exons'} = $exonCount;
	$sym2info{$geneSymbol}{'chr'} = $chrom;
	$sym2info{$geneSymbol}{'start'} = $txStart;
	$sym2info{$geneSymbol}{'end'} = $txEnd;

	my @es = split(",", $exonStarts);
	my @ee = split(",", $exonEnds);

	# transcript length
	{
		my $trlen;		
		for (my $i = 0; $i < @es; $i ++)
		{
			$trlen += $ee[$i]-$es[$i];
		}

		#$sym2info{$prev2sym{$geneSymbol}}{'trlen'} = $trlen if ($prev2sym{$geneSymbol});
		$sym2info{$geneSymbol}{'trlen'} = $trlen;
	}
	
	# compute cds length	
	if ($cdsStart and $cdsStart < $cdsEnd)
	{
		#print "$strand\t$cdsStart\t$cdsEnd\t$exonStarts\t$exonEnds\t";
		my ($st, $en, $cdslen);		
		for (my $i = 0; $i < @es and $cdsEnd > $es[$i]; $i ++)
		{
			next if ($cdsStart > $ee[$i]);
			$st = ($cdsStart > $es[$i] and $cdsStart < $ee[$i]) ? $cdsStart : $es[$i];
			$en = ($cdsEnd > $es[$i] and $cdsEnd < $ee[$i]) ? $cdsEnd : $ee[$i];
			$cdslen += $en-$st;
		}
	
		#$sym2info{$prev2sym{$geneSymbol}}{'cdslen'} = $cdslen if ($prev2sym{$geneSymbol});
		$sym2info{$geneSymbol}{'cdslen'} = $cdslen;
	}
}
close(G);
INFO("$lines genes read from file $ENV{HOME}/generic/data/hg19/hg19.knownGene.txt");

# TABLE: filtered-variants
my %genes;
my $written = 0;
<>; # skip header
while(<>)
{
	chomp;
	my ($patient, $sample, $var_type, $status, $rejected_because, $chr, $pos, $dbSNP, $ref, $alt, $gene, $add_genes, $impact, $effect, $non_silent, $deleterious, $exons, 
		$dp_rem_tot, $dp_rem_ref, $dp_rem_var, $freq_rem, $dp_leu_tot, $dp_leu_ref, $dp_leu_var, $freq_leu, $aa_change, $snpeff,
		$polyphen2, $sift, $gerp, $siphy, $interpro) = split("\t");

	die "ERROR: $0: snpeff annotation missing from following line:\n$_\n"
		if (!$snpeff);

	next if ($status eq "REJECT");
	next if ($effect =~ /^(DOWNSTREAM|INTERGENIC|INTRON|UPSTREAM|INTERGENIC_CONSERVED)$/);
	
	if ($effect !~ /^(SYNONYMOUS_START|SYNONYMOUS_CODING|SYNONYMOUS_STOP|UTR_5_PRIME|UTR_5_DELETED|START_GAINED|UTR_3_PRIME|UTR_3_DELETED|INTRON_CONSERVED|INTRAGENIC|EXON)$/)
	{
		$genes{$patient}{$sample}{$gene}{'nonsyn'}{"$chr:$pos:$ref->$alt:$freq_leu:$impact:$effect"} = $genes{$patient}{$sample}{$gene}{'nonsyn'}{"$chr:$pos:$ref->$alt:$freq_leu:$impact:$effect"} 
			? $genes{$patient}{$sample}{$gene}{'nonsyn'}{"$chr:$pos:$ref->$alt:$freq_leu:$impact:$effect"}.",$snpeff"
			: $snpeff;		

		$genes{$patient}{$sample}{$gene}{'max_af_ns'} = exists $genes{$patient}{$sample}{$gene}{'max_af_ns'} 
			? $genes{$patient}{$sample}{$gene}{'max_af_ns'} < $freq_leu ? $freq_leu : $genes{$patient}{$sample}{$gene}{'max_af_ns'}
			: $freq_leu;
		
		if ($exons)
		{
			my $exno = get_exon_no($exons, $gene);
			$genes{$patient}{$sample}{$gene}{'exons_ns'}{$exno} = 1;
		}
	}
	
	$genes{$patient}{$sample}{$gene}{'all'}{"$chr:$pos:$ref->$alt:$freq_leu:$impact:$effect"} = $genes{$patient}{$sample}{$gene}{'all'}{"$chr:$pos:$ref->$alt:$freq_leu:$impact:$effect"} 
		? $genes{$patient}{$sample}{$gene}{'all'}{"$chr:$pos:$ref->$alt:$freq_leu:$impact:$effect"}.",$snpeff"
		: $snpeff;

	$genes{$patient}{$sample}{$gene}{'max_af'} = exists $genes{$patient}{$sample}{$gene}{'max_af'} 
		? $genes{$patient}{$sample}{$gene}{'max_af'} < $freq_leu ? $freq_leu : $genes{$patient}{$sample}{$gene}{'max_af'}
		: $freq_leu;

	$genes{$patient}{$sample}{$gene}{'domains'} = defined $genes{$patient}{$sample}{$gene}{'domains'} 
		? $genes{$patient}{$sample}{$gene}{'domains'}."|$interpro"
		: $interpro;

	if ($exons)
	{
		my $exno = get_exon_no($exons, $gene);
		$genes{$patient}{$sample}{$gene}{'exons'}{$exno} = 1;
	}
}

# TABLE: impacted-genes
print "patient\tcomparison\tgene\tchr\tstart\tend\ttr_len\tcds_len\texons\tcosmic\tdesc\tnum_mut\tnum_mut_nonsyn\tmax_af\tmax_af_ns\timp_exons\timp_exons_ns\tmut_effects\tdomains\n";
foreach my $p (keys(%genes))
{
	foreach my $s (keys(%{$genes{$p}}))
	{
		my @sorted_genes = sort {values(%{$genes{$p}{$s}{$b}{'all'}}) <=> values(%{$genes{$p}{$s}{$a}{'all'}})} keys(%{$genes{$p}{$s}});
		foreach my $g (@sorted_genes)
		{
			my $info = $sym2info{$g}
				or WARN("Could not map gene $g");
						
			print $p, "\t", $s, "\t", $g, "\t";
			print "".($info->{'chr'} ? $info->{'chr'} : "")."\t";
			print "".($info->{'start'} ? $info->{'start'} : "")."\t";
			print "".($info->{'end'} ? $info->{'end'} : "")."\t";
			print "".($info->{'trlen'} ? $info->{'trlen'} : "")."\t";
			print "".($info->{'cdslen'} ? $info->{'cdslen'} : "")."\t";
			print "".($info->{'exons'} ? $info->{'exons'} : "")."\t";
			print $cosmic{$g} ? $cosmic{$g} : "", "\t"; 
			print "".($info->{'description'} ? $info->{'description'} : "")."\t";
			
			print scalar(values(%{$genes{$p}{$s}{$g}{'all'}})), "\t";
			print scalar(values(%{$genes{$p}{$s}{$g}{'nonsyn'}})), "\t";

			print $genes{$p}{$s}{$g}{'max_af'} ? $genes{$p}{$s}{$g}{'max_af'} : "", "\t";
			print $genes{$p}{$s}{$g}{'max_af_ns'} ? $genes{$p}{$s}{$g}{'max_af_ns'} : "", "\t";

			print $genes{$p}{$s}{$g}{'exons'} ? join(",", keys(%{$genes{$p}{$s}{$g}{'exons'}})) : "", "\t";
			print $genes{$p}{$s}{$g}{'exons_ns'} ? join(",", keys(%{$genes{$p}{$s}{$g}{'exons_ns'}})) : "", "\t";
			 
			my $first = 1;
			foreach my $v (keys(%{$genes{$p}{$s}{$g}{'all'}}))
			{
				print ";" if (!$first);
				print $v,"[",$genes{$p}{$s}{$g}{'all'}{$v},"]";
				$first = 0;
			}
			
			print "\t", $genes{$p}{$s}{$g}{'domains'} ? $genes{$p}{$s}{$g}{'domains'} : "";
			print "\n";
			
			$written ++;
		}		
	}
}
INFO("$written output lines written.");

my %chosen_transcript;
sub get_exon_no
{
	my $exons = shift or die "exon string not specified\n";
	my $gene = shift or die "gene not specified\n";

	my $exon_difftrans;	
	while($exons =~ /(\d+) \(([^,]+)\),?/g)
	{
		my ($exno, $genes) = ($1, $2);

		while ($genes =~ /([^:]+):([^;]+);?/g)
		{
			next if ($1 ne $gene);
#			print STDERR "EXNO: $exno, GENE: $1, TRANSCRIPT: $2\n" if ($gene eq "ATXN7");	
			if (exists $chosen_transcript{$gene} and $chosen_transcript{$gene} ne $2)
			{
				$exon_difftrans = $exno;
				next;
			}

			$chosen_transcript{$gene} = $2;
			return $exno;
		}
	}
	
	# this is actually a fallback if we could not find the chosen transcript for this gene
	return $exon_difftrans;
}