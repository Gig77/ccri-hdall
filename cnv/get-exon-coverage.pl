use warnings FATAL => qw( all );;
use strict;

use Getopt::Long;
use Set::IntervalTree;

my ($exon_bed_file);
GetOptions
(
	"exon-bed=s" => \$exon_bed_file
);

die "ERROR: --exon-bed not specified" if (!$exon_bed_file);

my %it;
open(B, "$exon_bed_file") or die "ERROR: Could not open file $exon_bed_file\n";
while(<B>)
{
	next if /^#/;
	next if /^$/;

	my ($chr, $start, $end) = split /\t/;

	#print "$chr,$start,$end\n";
	
	my $t = $it{$chr};
	if (!$t)
	{
		$t = Set::IntervalTree->new;
		$it{$chr} = $t;
	}

	my $cov = 0;
	$t->insert(\$cov, $start, $end);	
}
close(B);

while(<>)
{
        next if /^#/;
        next if /^$/;

	my ($chr, $pos, $cov) = split /\t/;

	my $t = $it{$chr};
	my $e = $t->fetch($pos, $pos+1) if ($t);
	if (!$e)
	{
		print STDERR "ERROR: Chromosomal coordinate not found in BED file: $chr:$pos\n";
		next;
	}

	${$e->[0]} += $cov;
}

open(B, "$exon_bed_file") or die "ERROR: Could not open file $exon_bed_file\n";
while(<B>)
{
        next if /^#/;
        next if /^$/;

	chomp;
        my ($chr, $start, $end) = split /\t/;

        my $t = $it{$chr};
        my $cov = $t->fetch($start, $end) if ($t);
        if (!$cov)
        {
                print STDERR "INTERNAL ERROR: BED region not found in interval tree: $_\n";
                next;
        }

	print $_, "\t", ${$cov->[0]}, "\t", ${$cov->[0]} / ($end-$start), "\n";
}       
close(B);
