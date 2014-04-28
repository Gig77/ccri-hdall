use warnings FATAL => qw( all );
use strict;

use Carp;
use Getopt::Long;

my ($format, $header);
GetOptions
(
	"format-dia=s" => \$format,
	"header" => \$header  
);

die "ERROR: Format not specified\n" if (!$format);

# TABLE: hdall.cnv.tsv
print "patient\tsample\tevent\tcnumber\tchromosome\tstart\tend\tsize\tnum_genes\tgenes\n" 
	if ($header);

<>; # skip header
if ($format eq 'andrea')
{
	while(<>)
	{
		chomp;
		my ($chromosome, $start, $end, $sample_id, $mean, $copy_number, $overlapping_features, $nearest_feature) = split /\t/;

		my ($patient, $sample) = $sample_id =~ /^(\d+)(D|R)$/;	
		if (!$patient or !$sample)
		{
			print STDERR "WARNING: Could not parse following line. SKIPPED.\n $_\n";
			next;
		}
		
		$sample =~ s/D/rem_dia/;
		$sample =~ s/R/rem_rel/;
		
		my %copynr2event = (
			'Amplification' => 'gain',
			'Deletion' => 'loss',
			'biDeletion' => 'loss',
			'UPD' => 'loh',
			'Tetrasom' => 'gain',
			'UPTri' => 'gain'
		);
		my $event = $copynr2event{$copy_number};
		die "ERROR: Could not map copy number to event: $copy_number\n" if (!$event);
		
		my $cnumber = 2;
		if ($mean < 0.6) { $cnumber = 0; }
		elsif ($mean < 1.6) { $cnumber = 1; }
		elsif ($mean >= 2.3 and $mean < 3.1) { $cnumber = 3; }
		elsif ($mean >= 3.1 and $mean < 4.6) { $cnumber = 4; }
		elsif ($mean >= 4.6 and $mean < 5.6) { $cnumber = 5; }
		elsif ($mean >= 5.6) { $cnumber = "6+"; }
		
		# correct glitch in input files
		$cnumber = 3 if ($event eq "gain" and $cnumber == 0);
		
		my $genes = $overlapping_features;
		$genes =~ s/region, ends, 0, bp, before, CNTNAP2, \(\+\)//;
		$genes =~ s/region overlaps with .+? of //g;
		$genes =~ s/region starts .+? after [^,]+//g;
		$genes =~ s/region ends .+? before [^,]+//g;
		$genes =~ s/contained within //g;
		$genes =~ s/intron of [^,]+//g;
		$genes =~ s/ (\(\-\)|\(\+\))//g;
		$genes =~ s/ , /, /g;

		print "$patient\t";
		print "$sample\t";
		print "$event\t";
		print "$cnumber\t";
		print "$chromosome\t";
		print "$start\t";
		print "$end\t";
		print "".($end-$start)."\t";
		print scalar(split(", ", $genes))."\t"; 
		print "$genes\n";
	}
}
elsif ($format eq 'maria')
{
	while(<>)
	{
		chomp;
		my ($file, $cn_state, $type, $chromosome, $min, $max, $size, $marker_count, $cytoband_start, $cytoband_end, $genes, $sno_mirna) = split /\t/;

		my ($patient, $sample) = $file =~ /^(C|B|\d{3})_?(Dx|DX|Rel|REL|RR|Rem)\./;
		if (!$patient or !$sample)
		{
			print STDERR "WARNING: Could not parse following line. SKIPPED.\n$_\n";
			next;
		}
		
		$sample =~ s/dx/rem_dia/i;
		$sample =~ s/rel/rem_rel/i;
		$sample =~ s/RR/rem_rel3/i; # patient 715

		
		if ($sample ne "rem_dia" and $sample ne "rem_rel" and $sample ne "rem_rel3")
		{
			print STDERR "WARNING: Uknown sample: $sample. SKIPPED.\n $_\n";
			next;
		}
		my $event = lc($type);
		$cn_state = 2 if (!defined $cn_state or $cn_state eq '');
		
		print "$patient\t";
		print "$sample\t";
		print "$event\t";
		print "$cn_state\t";
		print "$chromosome\t";
		print "$min\t";
		print "$max\t";
		print "".($max-$min)."\t";
		print scalar(split(", ", $genes))."\t"; 
		print "$genes\n";
	}
}
else
{
	die "ERROR: Invalid format: $format\n";
}