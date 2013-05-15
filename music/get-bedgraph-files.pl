#!/usr/bin/perl

use strict;
use LWP::Simple;

my $url = 'http://info.biomedical-sequencing.at/DNA/Aijae5aegah6ohph2yadau/bam/';
my $target_dir = "/data/christian/kamilla/data/current/bam/";

chdir($target_dir);

my $content = get($url);
while($content =~ /href=\"(.*?\.bedgraph)\"/g)
{
	my $cmd = "curl -O $url$1";
	print "$cmd\n";
	system("$cmd");
}
