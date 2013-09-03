use warnings FATAL => qw( all );
use strict;

use Carp;
use Getopt::Long;

# parse detailed results first
my ($sm_pathways, $sm_pathways_detail, $pathway_file, $maf_file);
GetOptions
(
	"pathway-file=s" => \$pathway_file, # input pathway file for music path-scan
	"sm-pathways=s" => \$sm_pathways, # result list from genome music path-scan  
	"sm-pathways-detail=s" => \$sm_pathways_detail, # detaild result from genome music path-scan (= second output file)  
	"maf-file=s" => \$maf_file # for determining type of mutation affecting each gene (needed to color KEGG diagrams)  
);

# read pathway input file
my %pathways;
open(P, "$pathway_file") or die "ERROR: could not open pathway file $pathway_file\n";
while(<P>)
{
	chomp;
	my ($path_id, $path_name, $class, $gene_line, $diseases, $drugs, $description) = split(/\t/);
	$pathways{$path_id}{'size'} = scalar(split('\|', $gene_line));
}
close(P);
croak "ERROR: File $pathway_file is empty" if (keys(%pathways) == 0);

my $content = "";
open(D,"$sm_pathways_detail") or die "ERROR: could not read file $sm_pathways_detail\n";
while(<D>) 
{
#	chomp;
	$content .= $_;
}
close(D);
croak "ERROR: File $sm_pathways_detail is empty" if ($content eq "");

# TABLE: cnv.maf
my %mut_type;
open(M,"$maf_file") or die "ERROR: could not read file $maf_file\n";
<M>;<M>; # skip headers
while(<M>) 
{
	chomp;
	my @f = split("\t");
	
	my $hugo = $f[0];
	my $center = $f[2];
	my $type = $f[9];
	
	if ($center eq 'Array')
	{
		if ($type eq "ins")
		{
			$mut_type{$hugo}{gain} = 1;
		}
		elsif ($type eq "del")
		{
			$mut_type{$hugo}{loss} = 1;
		}
		else
		{
			croak "ERROR: Invalid variant type $type:\n$_\n";
		}
	}
	else
	{
		$mut_type{$hugo}{mutated} = 1; 
	}
}
close(M);
croak "ERROR: File $maf_file is empty" if (keys(%mut_type) == 0);

my (%pw, %genes);
while ($content =~ /Pathway: (.*?)\nName: (.*?)\nClass: (.*?)\nP-value: (.*?)\nFDR: (.*?)\nDescription:(.*?)\n(.*?)\nSamples with mutations \(#hits\): (.*?)\n/sg)
{
#	next if ($1 ne 'GOTERM_BP_FAT:GO:0008284');
#	print "Pathway: $1\nName: $2\nClass: $3\nP-value: $4\nFDR: $5\nDescription: $6\n";
		
	my $id = $1;
	$pw{$id}{name} = $2;
	$pw{$id}{class} = $3;
	$pw{$id}{pvalue} = $4;
	$pw{$id}{fdr} = $5;
	$pw{$id}{desc} = $6;
	
	my $genes_per_sample = $7;
	foreach my $g (split(/\n/, $genes_per_sample))
	{
		my ($sample, $genes) = split (":", $g);
#		print "  $g\n";
		$pw{$id}{genes_per_sample}{$sample} = $genes;
		
		foreach my $gene (split(",", $genes))
		{
			$pw{$id}{genes}{$gene} = $pw{$id}{genes}{$gene} ? $pw{$id}{genes}{$gene} + 1 : 1;
			$genes{$gene}{$id} = 1;		
		}
	}
}	

# TABLE: sm_pathways.annotated
print "Pathway\tName\tClass\tSize\tSamples_Affected\tTotal_Variations\tp-value\tFDR\tNumGenes(deleted)\tGenes\tLink\tgenes_deleted\n";
open(D,"$sm_pathways") or die "ERROR: could not read file $sm_pathways\n";
<D>; # skip header
while(<D>) 
{
	chomp;
	my ($pathway, $name, $class, $num_samples, $num_genes, $pvalue, $fdr) = split(/\t/);
	
	my @genes_sorted = sort { $pw{$pathway}{genes}{$b} <=> $pw{$pathway}{genes}{$a} } keys(%{$pw{$pathway}{genes}});
	my @genes_print = map { "$_(".$pw{$pathway}{genes}{$_}.")"  } @genes_sorted; 

	# determine number of deleted genes
	my $num_deleted = 0;
	for (my $i = 0; $i < @genes_sorted; $i ++)
	{
		$num_deleted ++ if ($mut_type{$genes_sorted[$i]}{loss});
	}
	
	print "$pathway\t$name\t$class\t",$pathways{$pathway}{'size'},"\t$num_samples\t";
	print $num_genes, $num_deleted > 0 ? " ($num_deleted)" : "", "\t";
	print "$pvalue\t$fdr\t";
	print scalar(keys(%{$pw{$pathway}{genes}})),"\t";
	print join(",", @genes_print),"\t";
	if ($class  =~ /^KEGG/)
	{
		my $keggid = $pathway;
		$keggid =~ s/KEGG_PATHWAY://;
		print "http://www.genome.jp/kegg-bin/show_pathway?map=$keggid&multi_query=";
		for (my $i = 0; $i < @genes_sorted; $i ++)
		{
			my ($bgcolor, $fgcolor);
			if ($mut_type{$genes_sorted[$i]}{mutated} and $mut_type{$genes_sorted[$i]}{loss})
			{
				($bgcolor, $fgcolor) = ("yellow", "red");
			}
			elsif ($mut_type{$genes_sorted[$i]}{mutated} and $mut_type{$genes_sorted[$i]}{gain})
			{
				($bgcolor, $fgcolor) = ("yellow", "blue");
			}
			elsif ($mut_type{$genes_sorted[$i]}{loss})
			{
				($bgcolor, $fgcolor) = ("red", "black");
			}
			elsif ($mut_type{$genes_sorted[$i]}{gain})
			{
				($bgcolor, $fgcolor) = ("blue", "yellow");
			}
			else
			{
				($bgcolor, $fgcolor) = ("yellow", "black");
			}
			
			print $genes_sorted[$i],'+',$bgcolor,'%2C',$fgcolor;
			print '%0D%0A' if ($i < @genes_sorted - 1);
		}
	}
		
	print "\n";
}
close(D);
