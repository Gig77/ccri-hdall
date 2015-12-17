use warnings FATAL => qw( all );
use strict;
use Carp;

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
print STDERR "".(scalar(keys(%id2sym))." id mappgins read from file /mnt/projects/hdall/results/id-mappings.tsv\n");

# read biomart id mapping to get previous symbols
my (%ncbi2hugo, %ens2hugo);
open(G, "/mnt/projects/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv") or die "ERROR: could not read gene list\n";
while(<G>)
{
	chomp;
	my ($approved_symbol, $entrez_gene_id, $accession_numbers, $approved_name, 
		$previous_symbols, $previous_names, $aliases, $name_aliases, $entrez_gene_id_ncbi, $ensembl_id, $uniprot, $refseq, $ucsc) = split("\t");

	$approved_symbol = $id2sym{$approved_symbol} if ($approved_symbol and $id2sym{$approved_symbol});
	$ncbi2hugo{$entrez_gene_id} = $approved_symbol if ($entrez_gene_id and $approved_symbol);
	$ncbi2hugo{$entrez_gene_id_ncbi} = $approved_symbol if ($entrez_gene_id_ncbi and $approved_symbol);
	$ens2hugo{$ensembl_id} = $approved_symbol if ($ensembl_id and $approved_symbol);
}
close(G);

open(P, "/mnt/projects/hdall/results/music/pathscan/pathways-kegg.jerry.ensembl.tsv") or die "ERROR: could not read KEGG input file\n";
while(<P>)
{
	chomp;
	my ($pwid, $pwname, $pwcategory, $genes) = split("\t");
	$pwid =~ s/path://;
	print "$pwid\tKEGG:$pwname\tKEGG:$pwcategory\t";
	my @trgenes;
	foreach my $g (split('\|', $genes))
	{
		my ($ncbi, $ens) = split(":", $g);
		my $hugo = $ncbi2hugo{$ncbi};
		$hugo = $ens2hugo{$ens} if (!$hugo);
		$hugo = $ens if (!$hugo); # fallback
		push(@trgenes, "$ncbi:$hugo");
	}
	print join("\|", @trgenes);
	print "\t\t\t\n";
}
close(P);

