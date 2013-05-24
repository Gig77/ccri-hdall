#!/usr/bin/perl

use strict;
use warnings;

my $dir = "/home/STANNANET/christian.frech/hdall/data/bam";

chdir($dir);
opendir(D, $dir) or die "ERROR: could not read directory\n";
while (my $f = readdir(D))
{
	next if ($f !~ /bedgraph.gz$/);
	my $wigf = $f;
	$wigf =~ s/bedgraph.gz/music.wig/;
	next if (-s $wigf);
	
	next if ($wigf !~ /_rel/);
	my $cmd = "perl /home/STANNANET/christian.frech/hdall/scripts/bedgraph-to-wig.pl $f $wigf";
	print "$cmd\n";
	system("$cmd");
}
closedir(D);