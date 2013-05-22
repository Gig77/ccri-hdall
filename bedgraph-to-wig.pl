#!/usr/bin/perl
# script adapted from http://davetang.org/wiki/tiki-index.php?page=wig

use strict;
use warnings;

my $usage = "Usage: $0 <infile.bedgraph> <outfile.wig>\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;

if ($infile =~ /\.gz$/){
   open(IN,'-|',"gunzip -c $infile") || die "Could not open $infile: $!\n";
} else {
   open(IN,'<',$infile) || die "Could not open $infile: $!\n";
}

open(OUT,'>',$outfile) || die "Could not open $outfile: $!\n";
#print OUT "track type=wiggle_0 name=\"$infile\" description=\"$infile\" visibility=full\n";

while(<IN>){
   chomp;
   next if (/^track/);
   #chr1    3000403 3000404 2
   my ($chr,$start,$end,$data) = split(/\t/);
   $chr =~ s/^chr//;
   #fixedStep chrom=chr1 start=3016975 step=1 span=1
   #1
   #1
   #1
   my $length = $end - $start;
   print OUT "fixedStep chrom=$chr start=$start step=1 span=1\n";
   for (0 .. $length){
      print OUT "$data\n";
   }
#	print OUT "variableStep chrom=$chr span=",$end-$start,"\n";
#	print OUT "$start $data\n"; 
}
close(IN);
close(OUT);

exit(0);