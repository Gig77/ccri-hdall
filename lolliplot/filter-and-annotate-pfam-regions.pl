use strict;
use warnings FATAL => qw( all );

# read id mapping
my %id2sym;
open(M, "$ENV{HOME}/hdall/results/id-mappings.tsv") or die "ERROR: could not read id mappings\n";
while(<M>)
{
	chomp;
	my ($sym, $id) = split(/\t/);
	$id2sym{$id} = $sym;
}
close(M);

# uniprot to ucsc mappings
my %uniprot2refseq;
my %uniprot2ucsc;
my %uniprot2hugo;
open(G,"$ENV{HOME}/generic/data/hg19/hg19.kgXref.txt") or die "ERROR: could not open file $ENV{HOME}/generic/data/hg19/hg19.kgXref.txt";
while(<G>)
{
	chomp;
	my ($kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refSeq, $protAcc, $description, $rfamAcc, $tRnaName) = split(/\t/);

	$spID =~ s/\-\d+//;
	$uniprot2hugo{$spID} = $id2sym{$kgID};
	
	if ($mRNA =~ /^NM_/)
	{
		$uniprot2refseq{$spID} = $mRNA;
		$uniprot2ucsc{$spID} = $kgID;
	}
	elsif ($refSeq =~ /^NM_/)
	{
		$uniprot2refseq{$spID} = $refSeq;
		$uniprot2ucsc{$spID} = $kgID if (!$uniprot2ucsc{$spID});
	}	
}
close(G);

# read pfam domain descriptions
my %pfamid2name;
open(D, "$ENV{HOME}/generic/data/pfam-27.0/pfamA.txt") || die "ERROR: reading file $ENV{HOME}/generic/data/pfam-27.0/pfamA.txt";
while(<D>)
{
	chomp;
	my ($no, $pfamid, $name, $alias, $desc) = split("\t");
	$pfamid2name{$pfamid} = $name;
}
close(D);

open(D, "gunzip -c $ENV{HOME}/generic/data/pfam-27.0/Pfam-A.regions.tsv.gz |") || die "can't open pipe to $ENV{HOME}/generic/data/pfam-27.0/Pfam-A.regions.tsv.gz";
while(<D>)
{
	chomp;
	my ($id, $dummy1, $dummy2, $dummy3, $pfamid, $start, $end) = split("\t");

	next if (!defined $uniprot2hugo{$id} or !defined $uniprot2ucsc{$id} or !defined $uniprot2refseq{$id} or !defined $pfamid2name{$pfamid});	
	
	print $uniprot2hugo{$id}, "\t";
	print $uniprot2ucsc{$id}, "\t";
	print $uniprot2refseq{$id}, "\t$id\t$pfamid\t";
	print $pfamid2name{$pfamid}, "\t$start\t$end\n";
}
close(D);
