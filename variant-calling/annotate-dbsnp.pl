use warnings FATAL => qw( all );
use strict;

use Tabix;

my $dbsnp = Tabix->new(-data => $ENV{HOME}.'/generic/data/hg19/hg19.snp137Common.txt.gz');

while (<>)
{
	# skip header lines
	if (/^#/)
	{
		print $_;
		next;
	}
	
	my @f = split("\t");
	next if (length($f[4]) > 1 or length($f[5]) > 1); # ignore indels for now...
	my $iter = $dbsnp->query($f[0], $f[1]-1, $f[1]);

	if (!$iter or !$iter->{_})
	{
		print $_;
		next;
	}
	
	my @snps;
	while (my $line = $dbsnp->read($iter)) 
	{
		my @s = split("\t", $line);
		push(@snps, $s[4]) if ($s[11] eq 'single' and $s[2] == $f[1]-1 and $s[3] == $f[1]);
	}
	
	$f[2] = join(",", @snps);
	print join("\t", @f);	
}
