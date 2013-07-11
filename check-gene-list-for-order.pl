use warnings FATAL => qw( all );
use strict;

my %matrix;
open(B, "/home/christian/hdall/results/kamilla/final-list-gene-patient-matrix-v8.tsv") or die "error opening 2nd file: $!\n";
while(<B>)
{
	chomp;
	$matrix{$_} = 1 if ($_);
}
close(B);

my %studio;
open(A,"/home/christian/hdall/results/kamilla/final-list-design-studio.tsv") or die "error opening file\n";
<A>;
while(<A>)
{
	chomp;
	
	my ($target_region, $upstream_bases, $downstream_bases, $chr, $start, $end, $selected_targets, $total_targets, 
		$target_type, $selected_probes, $total_probes, $desired_probe_spacing, $coverage, $added, $labels, $design_warnings) = split /\t/;
	
	my $found;
	foreach my $g (keys(%matrix))
	{
		if ($labels =~ /^($g)(\d+)/)
		{
#			print "$1-$2\n";
			$found = $1;
		}
	}

	if ($found)
	{
		$studio{$found} = 1;
	}
	else
	{
		print "not in matrix: $labels\n";
		$studio{$labels} = 1;
	}
	
	
}
close(A);

foreach my $m (keys(%matrix))
{
	print "not in design studio: $m\n" if (!$studio{$m})
}