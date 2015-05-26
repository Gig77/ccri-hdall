
open(P, "/mnt/projects/hdall/results/panel-genes.tsv") or die "could not read panel genes\n";
my %genes;
while(<P>)
{
	chomp;
	$genes{$_} = 1;
}
close(P);

open(I, "/mnt/projects/hdall/results/kamilla/final-list-design-studio.tsv") or die "could not read design studio file\n";
while(<I>)
{
	chomp;
	my ($target_region, $upstream_bases, $downstream_bases, $chromosome, $start, $stop, $selected_targets, $total_targets, 
		$target_type, $selected_probes, $total_probes, $desired_probe_spacing, $coverage, $added, $labels, $design_warnings) = split /\t/;
		
	foreach my $g (keys(%genes))
	{
		if ($labels =~ /($g)(\d+)/)
		{
			$start -= 2; # add padding left for splice site mutation
			$stop += 2; # add padding left for splice site mutation
			print "$chromosome\t$start\t$stop\t$1\t$2\n";
		}
	}
}
close(I);
