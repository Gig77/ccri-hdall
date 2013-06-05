use strict;
use warnings FATAL => qw( all );

use lib '/mnt/suse/home/STANNANET/christian.frech/git/gms-core/lib/perl';

use Genome;
use Genome::Model::Tools::Graph::MutationDiagram;

my $md = new MutationDiagram(-annotation => 'bla');
