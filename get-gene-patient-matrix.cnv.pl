use warnings FATAL => qw( all );
use strict;

use lib "/mnt/projects/generic/scripts";
use Generic;
use Log::Log4perl qw(:easy);
use List::Util qw(min max);

# STDIN: list with genes copy-number variant genes in each patient
# STDOUT: gene/patient matrix indicating the copy number for each gene and patient

use Getopt::Long;

# parse detailed results first
my ($max_genes);
GetOptions
(
	"max-genes=i" => \$max_genes
);

die "ERROR: missing parameter --max-genes\n"
	if (!defined $max_genes);

# read copy-number info
my (%gene_info, %cnvs, %patients, %all_genes);
my $cnv_read = 0;
<>; # skip header: patient\tsample\tgene\tchr\tstart\tend\ttrlen\tcdslen\texons\tcosmic\tdescr\tcnumber\tevent\tevent_coordinate\tevent_size\tnum_genes\n
while(<>)
{
	chomp;
	my ($patient, $sample, $gene, $chr, $start, $end, $trlen, $cdslen, $exons, $cosmic, $descr, $cnumber, $event, $event_coordinate, $event_size, $num_genes) = split(/\t/);

	next if ($num_genes > $max_genes);
	
	$patients{$patient} = 1;
	$all_genes{$gene} = 1;

	$gene_info{'chr'}{$gene} = $chr;
	$gene_info{'start'}{$gene} = $start;
	$gene_info{'end'}{$gene} = $end;	
	$gene_info{'trlen'}{$gene} = $trlen;
	$gene_info{'cdslen'}{$gene} = $cdslen;
	$gene_info{'exons'}{$gene} = $exons;
	$gene_info{'cosmic'}{$gene} = $cosmic;
	$gene_info{'descr'}{$gene} = $descr;

	$cnvs{'cnumber'}{$sample}{$patient}{$gene} = $cnumber;
	$cnvs{'event'}{$sample}{$patient}{$gene} = $event;
	$cnvs{'event_coordinate'}{$sample}{$patient}{$gene} = $event_coordinate;
	
	$cnvs{'num_gain'}{$sample}{$gene} = defined $cnvs{'num_gain'}{$sample}{$gene} ? $cnvs{'num_gain'}{$sample}{$gene} + 1 : 1 if ($event eq 'gain'); 
	$cnvs{'num_loss'}{$sample}{$gene} = defined $cnvs{'num_loss'}{$sample}{$gene} ? $cnvs{'num_loss'}{$sample}{$gene} + 1 : 1 if ($event eq 'loss'); 
	$cnvs{'num_loh'}{$sample}{$gene} = defined $cnvs{'num_loh'}{$sample}{$gene} ? $cnvs{'num_loh'}{$sample}{$gene} + 1 : 1 if ($event eq 'loh'); 
	$cnvs{'num_events'}{$sample}{$gene} = defined $cnvs{'num_events'}{$sample}{$gene} ? $cnvs{'num_events'}{$sample}{$gene} + 1 : 1; 

	$cnv_read ++;
}
INFO("$cnv_read CNVs read.");

print "gene\tdescr\tchr\tstart\tend\texons\ttr_len\tcds_len\tcosmic\t";
print "freq-dia\tfreq-dia-gain\tfreq-dia-loss\tfreq-dia-loh";
map { print "\t$_-dia" } keys(%patients);
print "\tfreq-rel\tfreq-rel-gain\tfreq-rel-loss\tfreq-rel-loh";
map { print "\t$_-rel" } keys(%patients);
print "\n";

my @sorted_genes = sort { ($cnvs{'num_events'}{'rem_rel'}{$b} ? $cnvs{'num_events'}{'rem_rel'}{$b} : 0) <=> ($cnvs{'num_events'}{'rem_rel'}{$a} ? $cnvs{'num_events'}{'rem_rel'}{$a} : 0) } keys(%all_genes);

foreach my $g (@sorted_genes)
{
	print "$g\t";
	print $gene_info{'descr'}{$g},"\t";
	print $gene_info{'chr'}{$g},"\t";
	print $gene_info{'start'}{$g},"\t";
	print $gene_info{'end'}{$g},"\t";
	print $gene_info{'exons'}{$g},"\t",$gene_info{'trlen'}{$g},"\t",$gene_info{'cdslen'}{$g},"\t",$gene_info{'cosmic'}{$g},"\t";

	print $cnvs{'num_events'}{'rem_dia'}{$g} ? $cnvs{'num_events'}{'rem_dia'}{$g} : "0", "\t";
	print $cnvs{'num_gain'}{'rem_dia'}{$g} ? $cnvs{'num_gain'}{'rem_dia'}{$g} : "0", "\t";
	print $cnvs{'num_loss'}{'rem_dia'}{$g} ? $cnvs{'num_loss'}{'rem_dia'}{$g} : "0", "\t";
	print $cnvs{'num_loh'}{'rem_dia'}{$g} ? $cnvs{'num_loh'}{'rem_dia'}{$g} : "0", "\t";
	
	foreach my $p (keys(%patients))
	{
		if (defined $cnvs{'cnumber'}{'rem_dia'}{$p}{$g})
		{
			print $cnvs{'event'}{'rem_dia'}{$p}{$g};	
			print '|';
			print $cnvs{'cnumber'}{'rem_dia'}{$p}{$g};	
			print '|';
			print $cnvs{'event_coordinate'}{'rem_dia'}{$p}{$g};	
		}
		else
		{
			print " ";
		}
		print "\t";
	}

	print $cnvs{'num_events'}{'rem_rel'}{$g} ? $cnvs{'num_events'}{'rem_rel'}{$g} : "0", "\t";
	print $cnvs{'num_gain'}{'rem_rel'}{$g} ? $cnvs{'num_gain'}{'rem_rel'}{$g} : "0", "\t";
	print $cnvs{'num_loss'}{'rem_rel'}{$g} ? $cnvs{'num_loss'}{'rem_rel'}{$g} : "0", "\t";
	print $cnvs{'num_loh'}{'rem_rel'}{$g} ? $cnvs{'num_loh'}{'rem_rel'}{$g} : "0", "\t";
	
	foreach my $p (keys(%patients))
	{
		if (defined $cnvs{'cnumber'}{'rem_rel'}{$p}{$g})
		{
			print $cnvs{'event'}{'rem_rel'}{$p}{$g};	
			print ':';
			print $cnvs{'cnumber'}{'rem_rel'}{$p}{$g};	
			print '|';
			print $cnvs{'event_coordinate'}{'rem_rel'}{$p}{$g};	
		}
		else
		{
			print " ";
		}
		print "\t";
	}

	print "\n";
}
