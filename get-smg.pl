use warnings;
use strict;

my $sample = $ARGV[0] or die "ERROR: sample type not specified ('dia' or 'rel')\n";
die "ERROR: invalid sample type specified\n" if ($sample ne 'dia' and $sample ne 'rel');

my $min_freq = $ARGV[1] or die "ERROR: minimum gene mutation frequency across samples not specified\n";

# read gene-patient matrix
<STDIN>; # skip header
while(<STDIN>)
{
	chomp;
	my @c = split("\t");
	if ($sample eq 'dia' and $c[1] and $c[1] >= $min_freq)
	{
		print $c[0],"\n";
	}
	elsif ($sample eq 'rel' and $c[24] and $c[24] >= $min_freq)
	{
		print $c[0],"\n";	
	}
}