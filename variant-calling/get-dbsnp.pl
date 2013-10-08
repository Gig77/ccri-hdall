use warnings FATAL => qw( all );
use strict;

use Tabix;
use Getopt::Long;

my ($dbsnp, $region);
GetOptions
(
	"dbsnp=s" => \$dbsnp,  # Tabix indexed dbSNP file from UCSC (snp137Common)
	"region=s" => \$region # e.g. chr1:50000
);

my $snpdb = Tabix->new(-data => "$dbsnp");

my ($chr, $start, $end) = $region =~ /([^:]+):(\d+)\-?(\d*/;
$end=$start if (!defined $end);
my $snps = get_snps($chr, $start, $end);

foreach my $k (keys(%{$snps}))
{
	print "$k\n";
}


sub get_rois
{
	my ($chr, $start, $end) = @_;
	croak "ERROR: bad input parameters" if (!$chr or !$start or !$end);

	my %snps;

	my $iter = $roi->query($chr, $start, $end);
	return \%rois if (!$iter or !$iter->{_});

	while (my $line = $roi->read($iter)) 
	{
		my $snp = (split("\t", $line))[4];  
		$rois{$roi} = 1; 
	}

	return \%rois;
}
