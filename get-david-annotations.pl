#!perl
use strict;
use warnings FATAL => qw( all );

use SOAP::Lite;
use HTTP::Cookies;
use Data::Dumper;
use Carp;

# read gene symbols from biomart id mapping table
my (%name2id, %id2name);
open(GENES, "/home/STANNANET/christian.frech/generic/data/ensembl/gene-id-mapping.biomart-0.7.tsv") or die "ERROR: could not read gene list\n";
<GENES>;
while(<GENES>)
{
	chomp;
	my ($approved_symbol, $entrez_gene_id, $accession_numbers, $approved_name, 
		$previous_symbols, $previous_names, $aliases, $name_aliases, $entrez_gene_id_ncbi) = split("\t");

	next if (!$approved_symbol or $approved_symbol =~ /withdrawn/);

	$entrez_gene_id = $entrez_gene_id_ncbi if (!$entrez_gene_id);
	next if (!$entrez_gene_id);

	$name2id{$approved_symbol} = $entrez_gene_id;
	map {$name2id{$_} = $entrez_gene_id} split(", ", $previous_symbols);
	$id2name{$entrez_gene_id} = $approved_symbol;		
}
close(GENES);

my $soap = SOAP::Lite                             
	-> uri('http://service.session.sample')                
	-> proxy('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService',
			cookie_jar => HTTP::Cookies->new(ignore_discard=>1));

 #user authentication by email address
 #For new user registration, go to http://david.abcc.ncifcrf.gov/webservice/register.htm
 my $check = $soap->authenticate('christian.frech@ccri.at')->result;
  	print STDERR "User authentication: $check\n";

 die "ERROR: authentication failed\n"
 	if (lc($check) ne "true"); 


 #list all annotation category names
 # BBID,BIND,BIOCARTA,BLOCKS,CGAP_EST_QUARTILE,CGAP_SAGE_QUARTILE,CHROMOSOME,COG_NAME,COG_ONTOLOGY,CYTOBAND,DIP,EC_NUMBER,ENSEMBL_GENE_ID,ENTREZ_GENE_ID,
 # ENTREZ_GENE_SUMMARY,GENETIC_ASSOCIATION_DB_DISEASE,GENERIF_SUMMARY,GNF_U133A_QUARTILE,GENETIC_ASSOCIATION_DB_DISEASE_CLASS,GOTERM_BP_2,GOTERM_BP_1,
 # GOTERM_BP_4,GOTERM_BP_3,GOTERM_BP_FAT,GOTERM_BP_5,GOTERM_CC_1,GOTERM_BP_ALL,GOTERM_CC_3,GOTERM_CC_2,GOTERM_CC_5,GOTERM_CC_4,GOTERM_MF_1,GOTERM_MF_2,
 # GOTERM_CC_FAT,GOTERM_CC_ALL,GOTERM_MF_5,GOTERM_MF_FAT,GOTERM_MF_3,GOTERM_MF_4,HIV_INTERACTION_CATEGORY,HIV_INTERACTION_PUBMED_ID,GOTERM_MF_ALL,
 # HIV_INTERACTION,KEGG_PATHWAY,HOMOLOGOUS_GENE,INTERPRO,OFFICIAL_GENE_SYMBOL,NCICB_CAPATHWAY_INTERACTION,MINT,PANTHER_MF_ALL,PANTHER_FAMILY,
 # PANTHER_BP_ALL,OMIM_DISEASE,PFAM,PANTHER_SUBFAMILY,PANTHER_PATHWAY,PIR_SUPERFAMILY,PIR_SUMMARY,PIR_SEQ_FEATURE,PROSITE,PUBMED_ID,REACTOME_INTERACTION,
 # REACTOME_PATHWAY,PIR_TISSUE_SPECIFICITY,PRINTS,PRODOM,PROFILE,SMART,SP_COMMENT,SP_COMMENT_TYPE,SP_PIR_KEYWORDS,SCOP_CLASS,SCOP_FAMILY,
 # SCOP_FOLD,SCOP_SUPERFAMILY,UP_SEQ_FEATURE,UNIGENE_EST_QUARTILE,ZFIN_ANATOMY,UP_TISSUE,TIGRFAMS,SSF,UCSC_TFBS
# my $allCategoryNames= $soap ->getAllAnnotationCategoryNames()->result;	 	  	
# print STDERR  "\nAll available annotation category names: \n$allCategoryNames\n";
 
 my $default_categories = "KEGG_PATHWAY,OMIM_DISEASE,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,BBID,BIOCARTA,INTERPRO,PIR_SUPERFAMILY,SMART";
 my $additional_categories = "UP_SEQ_FEATURE";
 $soap->setCategories("$default_categories,$additional_categories")
 	or croak "ERROR: could not set categories\n";
# print "$validated_categories[0]\n";
# exit;
 
 #addList
 #my $inputIds = 'KRAS,WSCD1,KRT17,MLL3,CREBBP,MGC70870,PDE4DIP,DEGS2,MST1L,DNAH9,XG,HOXB9,ZNF492,PDHA1';
 
 # get pathway annotation from DAVID
 # ATTENTION: we have to query in batches, otherwise we get into troubles with their sysadmins ...
 my %terms;
 my @ids = keys(%id2name);
 #@ids = reverse sort(@ids); # comment this out on second run! 
 while (my @sids = splice(@ids, 0, 500))
 {
	 my $inputIds = join(",", @sids);
	 #my $inputIds = 8563; # geht net
	 #my $inputIds = 57573; # geht
	 #my $inputIds = "8563,57573"; # get net
	 my $idType = 'ENTREZ_GENE_ID';
	 my $listName = 'cancer';
	 my $listType=0;
	 
	 print STDERR "Input genes: $inputIds", "\n";
	 #to add background list, set listType=1
	 my $list = $soap->addList($inputIds, $idType, $listName, $listType)->result;
	 print STDERR "Percentage of mapped genes: $list\n"; 
	 
	 #list all species  names
	# my $allSpecies= $soap ->getSpecies()->result;	 	  	
	 # print STDERR  "\nAll species: \n$allSpecies\n"; 
	 #list current species  names
	# my $currentSpecies= $soap->getCurrentSpecies()->result;	 	  	
	# print STDERR  "Current species: $currentSpecies\n"; 
	
	 #set user defined species 
	 #my $species = $soap ->setCurrentSpecies("1")->result;
	
	 #print STDERR "\nCurrent species: \n$species\n"; 
	
	#my $tableReport = $soap->getListReport();
	#print Dumper($tableReport->paramsout),"\n";
	#exit; 
	
	# MuSiC required format:
	# This is a tab-delimited file prepared from a pathway database (such as KEGG), with the columns: 
	# [path_id, path_name, class, gene_line, diseases, drugs, description] The latter three columns 
	# are optional (but are available on KEGG). The gene_line contains the "entrez_id:gene_name" of 
	# all genes involved in this pathway, each separated by a "|" symbol.
	# Example:
	# hsa00061      Fatty acid biosynthesis	      Lipid Metabolism     31:ACACA|32:ACACB|27349:MCAT|2194:FASN|54995:OXSM|55301:OLAH
	
	my $tableReport = $soap->getTableReport();
	#print Dumper($tableReport->paramsout),"\n";
	my $num_entries = 0;
	foreach my $r ($tableReport->paramsout)
	{
		my $entrez_id = ref($r->{'values'}) eq "ARRAY" ? $r->{'values'}->[0]->{'array'} : $r->{'values'}->{'array'};

		my @a_array = ref($r->{'annotationRecords'}) eq "ARRAY" ? @{$r->{'annotationRecords'}} : ($r->{'annotationRecords'});	
		foreach my $a (@a_array)
		{
			my $category = $a->{'category'};
			
			my @term_array = ref($a->{'terms'}) eq "ARRAY" ? @{$a->{'terms'}} : ($a->{'terms'});
			foreach my $t (@term_array)
			{
				my ($term_id, $term_name) = $t =~ /([^\$]+)\$(.+)/;
				if ($term_name =~ /(GO:[^\~]+)\~(.+)/)
				{
					($term_id, $term_name) = ($1, $2); 
				}
				elsif ($term_name =~ /([^:]+):(.+)/)
				{
					($term_id, $term_name) = ($1, $2); 
				}			
	
	#			print "$entrez_id\t$category\t$term_id\t$term_name\n";
							
				$terms{$term_id}{'name'} = $term_name;
				$terms{$term_id}{'category'} = $category;
				$terms{$term_id}{'genes'}{$entrez_id} = 1;
			}
		}
		$num_entries ++;	
	}
	print STDERR "Number of annotated entries: $num_entries\n";	
}

# output pathways and genes contained therein
foreach my $t (keys(%terms))
{
	print $terms{$t}{'category'},":$t\t", $terms{$t}{'name'},"\t",$terms{$t}{'category'},"\t";
	my @genes;
	map { push(@genes, "$_:".$id2name{$_}) } keys(%{$terms{$t}{'genes'}});
	print join("\|", @genes),"\t\t\t\n";
}

