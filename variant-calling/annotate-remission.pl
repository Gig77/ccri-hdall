use warnings FATAL => qw( all );
use strict;

use Tabix;

use Getopt::Long;

my ($rem_sample, $patient);
GetOptions
(
	"rem-sample=s" => \$rem_sample, # tabix indexed remission sample
	"patient=s" => \$patient  
);

die "--rem-sample not specified\n" if (!$rem_sample);
die "--patient not specified\n" if (!$patient);

my ($name_rem, $name_rel) = ("$patient"."_rem", "$patient"."_rel");
my $dbsnp = Tabix->new(-data => $rem_sample);

my ($total, $num_rem) = (0, 0);
while (<>)
{
	chomp;
	
	# skip header lines
	if (/^#/)
	{
		if (/^#CHROM/)
		{
			s/Sample1/$name_rel\t$name_rem/;
		}
		
		print $_, "\n";
		
		print '##FILTER=<ID=REM,Description="Variant is present in remission sample(s)">', "\n"
			if (/^##FILTER=<ID=indelError/);
		
		next;
	}
	
	$total ++;
	
	my @f = split("\t");

	my $iter = $dbsnp->query($f[0], $f[1]-1, $f[1]);
	if (!$iter or !$iter->{_})
	{
		print $_, "\n";
		next;
	}
	
	my $fr_line = $dbsnp->read($iter);
	if (!$fr_line)
	{
		print $_, "\n";
		next;
	}
	
	my @fr = split("\t", $fr_line);
	if ($fr[0] ne $f[0] or $fr[1] ne $f[1] or $fr[3] ne $f[3] or $fr[4] ne $f[4])
	{
		print $_, "\n";
		next;		
	}

	##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
	##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Raw Read Depth as reported by SAMtools">
	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= 15">
	##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
	##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
	##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
	##FORMAT=<ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
	##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description="Average quality of reference-supporting bases (qual1)">
	##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average quality of variant-supporting bases (qual2)">
	##FORMAT=<ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">
	##FORMAT=<ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">
	##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Depth of variant-supporting bases on forward strand (reads2plus)">
	##FORMAT=<ID=ADR,Number=1,Type=Integer,Description="Depth of variant-supporting bases on reverse strand (reads2minus)">

	my ($gt_t, $gq_t, $sdp_t, $dp_t, $rd_t, $ad_t, $freq_t, $pval_t, $rbq_t, $abq_t, $rdf_t, $rdr_t, $adf_t, $adr_t) = split(":", $f[9]);
	my ($gt_r, $gq_r, $sdp_r, $dp_r, $rd_r, $ad_r, $freq_r, $pval_r, $rbq_r, $abq_r, $rdf_r, $rdr_r, $adf_r, $adr_r) = split(":", $fr[9]);
	$freq_t =~ s/\%//;
	$freq_r =~ s/\%//;
	
	$f[6] = "REM"; # if ($freq_t < 30);
	$num_rem ++;
	
	print join("\t", @f), "\t", $fr[9], "\n";
}

print STDERR "$num_rem of $total variants present in remission sample(s)\n";