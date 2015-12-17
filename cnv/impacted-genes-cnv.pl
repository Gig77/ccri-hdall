use warnings FATAL => qw( all );
use strict;

use lib "/mnt/projects/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use Carp;
use Getopt::Long;

my $debug = 1;
$| = 1; # turn on autoflush

my ($max_genes);
GetOptions
(
	"max-genes=i" => \$max_genes
);
die("ERROR: --max-genes not specified\n") if (!defined $max_genes);

# read id mapping
my %id2sym;
open(M, "/mnt/projects/hdall/results/id-mappings.tsv") or croak "ERROR: could not read id mappings\n";
while(<M>)
{
	chomp;
	my ($sym, $id) = split(/\t/);
	$id2sym{$id} = $sym;
}
close(M);
INFO(scalar(keys(%id2sym))." id mappgins read from file /mnt/projects/hdall/results/id-mappings.tsv");

sub getsym 
{ 
	my $sym = shift;
	my $depth = 0; 
	while($id2sym{$sym} and $depth++<2) { $sym = $id2sym{$sym}}; 
	return $sym; 
}

# read gene description
my %sym2info;
open(G,"/mnt/projects/generic/data/hg19/hg19.kgXref.txt") or die "could not open file /mnt/projects/generic/data/hg19/hg19.kgXref.txt";
while(<G>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refSeq, $protAcc, $description, $rfamAcc, $tRnaName) = split(/\t/);

	#my $sym = $id2sym{$kgID};
	#$sym = $id2sym{$sym} if ($id2sym{$sym});
	my $sym = getsym($kgID);
	$sym2info{$sym}{'description'} = $description;	
}
close(G);
INFO(scalar(keys(%sym2info))." gene descriptions read from file /mnt/projects/generic/data/hg19/hg19.kgXref.txt");

# read cosmic
my $cosmic_file = "/mnt/projects/generic/data/cosmic/v65/cancer_gene_census.tsv"; 
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
open(G,"/mnt/projects/generic/data/hg19/hg19.knownCanonical.txt") or die "could not open file /mnt/projects/generic/data/hg19/hg19.knownCanonical.txt";
while(<G>)
{
	chomp;
	my ($chrom, $chromStart, $chromEnd, $clusterId, $transcript, $protein) = split(/\t/);
	my $geneSymbol = getsym($transcript);
	
	my $size = $chromEnd-$chromStart;
	
	# if multiple canonical transcripts for this gene symbol, use larger one
	next if ($canonical{$geneSymbol} and $sym2size{$geneSymbol} > $size); 
	
	$canonical{$geneSymbol} = $transcript;
	$sym2size{$geneSymbol} = $size;
}
close(G);
INFO(scalar(keys(%canonical))." canonical genes read from file /mnt/projects/generic/data/hg19/hg19.knownCanonical.txt");

my $lines = 0;
open(G,"/mnt/projects/generic/data/hg19/hg19.knownGene.txt") or die "could not open file /mnt/projects/generic/data/hg19/hg19.knownGene.txt";
while(<G>)
{
	chomp;
	my ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $proteinID, $alignID) = split(/\t/);

	$lines++;
	my $geneSymbol = getsym($name);
	
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
INFO("$lines genes read from file /mnt/projects/generic/data/hg19/hg19.knownGene.txt");

my %imp_genes;
my $written = 0;
<>; # skip header
while(<>)
{
	chomp;
	my ($patient, $sample, $event, $cnumber, $chromosome, $start, $end, $size, $num_genes, $genes) = split("\t");

	die "ERROR: Could not parse following line:\n$_\n"
		if (!defined $genes);

	if ($num_genes > $max_genes)
	{
		INFO("$patient: $sample: $chromosome:$start..$end: More than $max_genes genes. Skipped.\n");
		next;
	}
	
	$genes =~ s/(\S),(\S)/$1, $2/g;
	foreach my $g (split(", ", $genes))
	{
		$g =~ s/\s+$//;
		$g =~ s/^\s+//;
		my $g = getsym($g);
		
		$imp_genes{$patient}{$sample}{$g}{'cnumber'} = $cnumber;
		$imp_genes{$patient}{$sample}{$g}{'event'} = $event;
		$imp_genes{$patient}{$sample}{$g}{'event_coord'} = "$chromosome:$start-$end";
		$imp_genes{$patient}{$sample}{$g}{'event_size'} = $end-$start;
		$imp_genes{$patient}{$sample}{$g}{'num_genes'} = $num_genes;
	}
}	

my %not_map;
print "patient\tsample\tgene\tchr\tstart\tend\ttrlen\tcdslen\texons\tcosmic\tdescr\tcnumber\tevent\tevent_coordinate\tevent_size\tnum_genes\n";
foreach my $p (keys(%imp_genes))
{
	foreach my $s (keys(%{$imp_genes{$p}}))
	{
		foreach my $g (keys(%{$imp_genes{$p}{$s}}))
		{		
			next if ($g eq "");	
			my $info = $sym2info{$g}
				or $not_map{$g} = 1;
						
			print $p, "\t", $s, "\t", $g, "\t";
			print "".($info->{'chr'} ? $info->{'chr'} : "")."\t";
			print "".($info->{'start'} ? $info->{'start'} : "")."\t";
			print "".($info->{'end'} ? $info->{'end'} : "")."\t";
			print "".($info->{'trlen'} ? $info->{'trlen'} : "")."\t";
			print "".($info->{'cdslen'} ? $info->{'cdslen'} : "")."\t";
			print "".($info->{'exons'} ? $info->{'exons'} : "")."\t";
			print $cosmic{$g} ? $cosmic{$g} : "", "\t"; 
			print "".($info->{'description'} ? $info->{'description'} : "")."\t";
			
			print $imp_genes{$p}{$s}{$g}{'cnumber'}, "\t";
			print $imp_genes{$p}{$s}{$g}{'event'}, "\t";
			print $imp_genes{$p}{$s}{$g}{'event_coord'}, "\t";
			print $imp_genes{$p}{$s}{$g}{'event_size'}, "\t";
			print $imp_genes{$p}{$s}{$g}{'num_genes'}, "\n";
			
			$written ++;
		}		
	}
}

WARN("Could not map following genes: ", join(", ", keys(%not_map)))
	if (keys(%not_map) > 0);
	
INFO("$written output lines written.");
