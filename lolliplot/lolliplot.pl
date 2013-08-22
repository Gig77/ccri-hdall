use strict;
use warnings FATAL => qw( all );

use lib '/home/christian/hdall/scripts/lolliplot';
use lib "$ENV{HOME}/generic/scripts";

use Carp;
use Generic;
use Log::Log4perl qw(:easy);
use Lolliplot;

my $hugos = "CREBBP,KRAS,NRAS,TRRAP,CDC42EP1,SCN5A,CTBS,FRG1,TBP";

# ucsc/HUGO mapping
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

# refseq/ucsc mapping
my %transcripts;
my %refseq2ucsc;
open(G,"$ENV{HOME}/hdall/data/hg19/hg19.kgXref.txt") or die "ERROR: could not open file $ENV{HOME}/hdall/data/hg19/hg19.kgXref.txt";
while(<G>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refSeq, $protAcc, $description, $rfamAcc, $tRnaName) = split(/\t/);

	if ($mRNA =~ /^NM_/)
	{
		$transcripts{$kgID}{'refseq'} = $mRNA;
		$refseq2ucsc{$mRNA} = $kgID;
	}
	elsif ($refSeq =~ /^NM_/)
	{
		$transcripts{$kgID}{'refseq'} = $refSeq;
		$refseq2ucsc{$refSeq} = $kgID;
	}	
}
close(G);

my $lines = 0;
open(G,"$ENV{HOME}/hdall/data/hg19/hg19.knownGene.txt") or die "could not open file $ENV{HOME}/hdall/data/hg19/hg19.knownGene.txt";
while(<G>)
{
	chomp;
	my ($ucsc_id, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $proteinID, $alignID) = split(/\t/);

	# compute protein length	
	my @es = split(",", $exonStarts);
	my @ee = split(",", $exonEnds);
	if ($cdsStart and $cdsStart < $cdsEnd)
	{
		my ($st, $en, $cdslen);		
		for (my $i = 0; $i < @es and $cdsEnd > $es[$i]; $i ++)
		{
			next if ($cdsStart > $ee[$i]);
			$st = ($cdsStart > $es[$i] and $cdsStart < $ee[$i]) ? $cdsStart : $es[$i];
			$en = ($cdsEnd > $es[$i] and $cdsEnd < $ee[$i]) ? $cdsEnd : $ee[$i];
			$cdslen += $en-$st;
			$transcripts{$ucsc_id}{'splice_pos'}{$i+1} = int($cdslen / 3);	
			print "SPLICE SITE ", $i+1, ": ", int($cdslen / 3), "\n" if ($ucsc_id eq "uc003xyq.3");	
		}
	
		$transcripts{$ucsc_id}{'protlen'} = $cdslen/3;
		$transcripts{$ucsc_id}{'hugo'} = $id2sym{$ucsc_id};
	}
	
	$lines++;
}
close(G);
INFO("$lines genes read from file $ENV{HOME}/hdall/data/hg19/hg19.knownGene.txt");

open(D,"$ENV{HOME}/hdall/results/lolliplot/pfam-regions.filtered.tsv") or die "could not open file $ENV{HOME}/hdall/results/lolliplot/pfam-regions.filtered.tsv";
my %domains;
while(<D>)
{
	chomp;
	my ($hugo, $ucsc, $refseq, $uniprot, $domain_id, $pfamname, $start, $end) = split(/\t/);
	
	$domains{"$ucsc:$domain_id"} = defined $domains{"$ucsc:$domain_id"} ? $domains{"$ucsc:$domain_id"} + 1 : 1;
	my $domain_key = "$domain_id:".$domains{"$ucsc:$domain_id"};
	
	$domains{$ucsc}{$domain_key}{id} = $domain_id;
	$domains{$ucsc}{$domain_key}{source} = "Pfam27";
	$domains{$ucsc}{$domain_key}{name} = $pfamname;
	$domains{$ucsc}{$domain_key}{start} = $start;
	$domains{$ucsc}{$domain_key}{end} = $end;
	
	$lines++;
}
close(D);
INFO("$lines domain annotations read from file $ENV{HOME}/hdall/results/lolliplot/pfam-regions.filtered.tsv");

# read variants
# TABLE: filtered-variants
$lines = 0;
my %variants;
open(V, "$ENV{HOME}/hdall/results/filtered-variants.cosmic.tsv") or die "could not open file $ENV{HOME}/hdall/results/filtered-variants.cosmic.tsv";
<V>; # skip header
while(<V>)
{
	chomp;
	my ($patient, $sample, $var_type, $status, $chr, $pos, $dbSNP, $ref, $alt, $gene, $add_genes, $impact, $effect_notused, $exons, 
		$dp_rem_tot, $dp_rem_ref, $dp_rem_var, $freq_rem, $dp_leu_tot, $dp_leu_ref, $dp_leu_var, $freq_leu, $aa_change, $snpeff,
		$polyphen2, $sift, $gerp, $siphy, $interpro, $af_1000g) = split("\t");

	$snpeff =~ s/EFF=//;
	foreach my $eff (split(",", $snpeff))
	{
		my ($effect, $rest) = $eff =~ /([^\(]+)\(([^\)]+)\)/
			or die "ERROR: could not parse SNP effect: $snpeff\n";

		# skip mutations not affecting coding sequence of protein			
		next if ($effect =~ /^(UPSTREAM|DOWNSTREAM|UTR_5_PRIME|UTR_5_DELETED|START_GAINED|UTR_3_PRIME|UTR_3_DELETED|INTRON|INTRON_CONSERVED|INTERGENIC|INTERGENIC_CONSERVED|INTRAGENIC|EXON)$/);

		my ($impact, $class, $codon_change, $aa_change, $aa_length, $gene_name, $gene_biotype, 
			$coding, $transcript, $exon, $genotype_num) = split('\|', $rest)
				or die "ERROR: could not parse SNP effect: $eff\n";

		$transcript =~ s/\.\d+$//; # remove version number from accession
		$transcript =~ s/\.\d+$//; 

		my $ucsc = $refseq2ucsc{$transcript};
		if (!$ucsc)
		{
			print STDERR "WARNING: Could not map RefSeq ID $transcript to UCSC ID.\n";
			next;
		}
		
		# map splice site mutations onto protein sequence
		if ($effect eq "SPLICE_SITE_ACCEPTOR")
		{
			$aa_change = $transcripts{$ucsc}{'splice_pos'}{$exon-1};
		}
		elsif ($effect eq "SPLICE_SITE_DONOR")
		{
			$aa_change = $transcripts{$ucsc}{'splice_pos'}{$exon};			
		}
		
		if (!$aa_change)
		{
			print "WARNING: Could not map following mutation to protein sequence:\n";
			print "$transcript $eff\n";
			next;
		}
		
		$variants{"dia"}{$ucsc}{"$patient\t$var_type\t$status\t$chr\t$pos\t$aa_change"} = "$effect" if ($sample eq "rem_dia");
		$variants{"rel"}{$ucsc}{"$patient\t$var_type\t$status\t$chr\t$pos\t$aa_change"} = "$effect" if ($sample eq "rem_rel");
		$variants{"both"}{$ucsc}{"$patient\t$var_type\t$status\t$chr\t$pos\t$aa_change"} = "$effect";		 
	}

	$lines ++;
}
close(V);
INFO("$lines variants read from file $ENV{HOME}/hdall/results/filtered-variants.cosmic.tsv");


Lolliplot->new
(
	hugos => $hugos,
	variants => $variants{"dia"}, 
	transcripts => \%transcripts, 
	domains => \%domains,
	output_directory => "$ENV{HOME}/hdall/results/lolliplot",
	basename => "lolliplot_",
	suffix => "_dia",
	lolli_shape => "circle"
);

Lolliplot->new
(
	hugos => $hugos,
	variants => $variants{"rel"}, 
	transcripts => \%transcripts, 
	domains => \%domains,
	output_directory => "$ENV{HOME}/hdall/results/lolliplot",
	basename => "lolliplot_",
	suffix => "_rel",
	lolli_shape => "circle"
);

Lolliplot->new
(
	hugos => $hugos,
	variants => $variants{"both"}, 
	transcripts => \%transcripts, 
	domains => \%domains,
	output_directory => "$ENV{HOME}/hdall/results/lolliplot",
	basename => "lolliplot_",
	suffix => "_both",
	lolli_shape => "circle"
);