#!/usr/bin/perl

use strict;
use LWP::Simple;

my $url = 'http://info.biomedical-sequencing.at/DNA/Aijae5aegah6ohph2yadau/bam/';
my $target_dir = "/home/STANNANET/christian.frech/hdall/data/bam/";

chdir($target_dir);

my $content = get($url);
while($content =~ /href=\"(.*?\.bedgraph)\"/g)
{
	if (!-e "$1" and !-e "$1.gz")
	{
		my $cmd = "curl -O $url$1";
		print "$cmd\n";		
		system("$cmd");
		system("gzip $1")
	}
}
