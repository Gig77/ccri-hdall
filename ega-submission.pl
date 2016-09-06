use warnings;
use strict;

use Data::Dumper;
use XML::Compile::Schema;
use XML::LibXML::Reader;

my $input_subjects = "/mnt/projects/hdall/results/ega-submission/input/subjects.csv";
my $input_samples = "/mnt/projects/hdall/results/ega-submission/input/samples.csv";
my $input_experiments = "/mnt/projects/hdall/results/ega-submission/input/experiments.csv";
my $input_analyses = "/mnt/projects/hdall/results/ega-submission/input/analyses.csv";
my $input_submission_result = "/mnt/projects/hdall/results/ega-submission/input/submission-result.txt";

my $output_sample_xml = "/mnt/projects/hdall/results/ega-submission/sample.xml";
my $output_experiment_xml = "/mnt/projects/hdall/results/ega-submission/experiment.xml";
my $output_run_xml = "/mnt/projects/hdall/results/ega-submission/run.xml";
my $output_analysis_xml = "/mnt/projects/hdall/results/ega-submission/analysis.xml";
my $output_dataset_xml = "/mnt/projects/hdall/results/ega-submission/dataset.xml";
my $output_submission_xml = "/mnt/projects/hdall/results/ega-submission/submission.xml";

my $schema_common = "ftp://ftp.sra.ebi.ac.uk/meta/xsd/sra_1_5/SRA.common.xsd";
my $schema_sample = "ftp://ftp.sra.ebi.ac.uk/meta/xsd/sra_1_5/SRA.sample.xsd";
my $schema_experiment = "ftp://ftp.sra.ebi.ac.uk/meta/xsd/sra_1_5/SRA.experiment.xsd";
my $schema_run = "ftp://ftp.sra.ebi.ac.uk/meta/xsd/sra_1_5/SRA.run.xsd";
my $schema_analysis = "ftp://ftp.sra.ebi.ac.uk/meta/xsd/sra_1_5/SRA.analysis.xsd";
my $schema_dataset = "ftp://ftp.sra.ebi.ac.uk/meta/xsd/sra_1_5/EGA.dataset.xsd";
my $schema_submission = "ftp://ftp.sra.ebi.ac.uk/meta/xsd/sra_1_5/SRA.submission.xsd";

#----------------
# download XML schemas
#----------------

#system("wget $schema_common -O /mnt/projects/hdall/results/ega-submission/schemas/SRA.common.xsd");
#system("wget $schema_sample -O /mnt/projects/hdall/results/ega-submission/schemas/SRA.sample.xsd");
#system("wget $schema_experiment -O /mnt/projects/hdall/results/ega-submission/schemas/SRA.experiment.xsd");
#system("wget $schema_run -O /mnt/projects/hdall/results/ega-submission/schemas/SRA.run.xsd");
#system("wget $schema_analysis -O /mnt/projects/hdall/results/ega-submission/schemas/SRA.analysis.xsd");
#system("wget $schema_dataset -O /mnt/projects/hdall/results/ega-submission/schemas/EGA.dataset.xsd");
#system("wget $schema_submission -O /mnt/projects/hdall/results/ega-submission/schemas/SRA.submission.xsd");

#----------------
# build data structure
#----------------

# read subject table

my %subjects;
open(SUBJECTS, "$input_subjects") or die "ERROR: Could not open file $input_subjects\n";
<SUBJECTS>; # skip header
while(<SUBJECTS>) {
	chomp;
	my ($subject_id, $tax_id, $scientific_name, $common_name, $cohort, $gender, $age_dia_years) = split("\t", $_);
	$subjects{$subject_id}{ID} = $subject_id;
	$subjects{$subject_id}{TAXON_ID} = $tax_id;
	$subjects{$subject_id}{SCIENTIFIC_NAME} = $scientific_name;
	$subjects{$subject_id}{COMMON_NAME} = $common_name;
	$subjects{$subject_id}{COHORT} = $cohort;
	$subjects{$subject_id}{GENDER} = $gender;
	$subjects{$subject_id}{AGE_AT_DIAGNOSIS_YEARS} = $age_dia_years;
}
close(SUBJECTS);

# read sample table

my %samples;
open(SAMPLES, "$input_samples") or die "ERROR: Could not open file $input_samples\n";
<SAMPLES>;  # skip header
while (<SAMPLES>) {
	chomp;
	my ($alias, $title, $description, $phenotype, $subject_id) = split("\t", $_);
	
	$samples{$alias}{ID} = $alias;
	$samples{$alias}{TITLE} = $title;
	$samples{$alias}{DESCRIPTION} = $description;
	$samples{$alias}{PHENOTYPE} = $phenotype;
	$samples{$alias}{SUBJECT} = $subjects{$subject_id};
	
	die "ERROR: Subject ID '$subject_id' not found in file $input_subjects\n" if (!exists $subjects{$subject_id});
	
	$subjects{$subject_id}{SAMPLES} = {} if (!exists $subjects{$subject_id}{SAMPLES});
	$subjects{$subject_id}{SAMPLES}{$alias} = $samples{$alias};
}
close(SAMPLES);

# read experiments

my %experiments;
open(EXPERIMENTS, "$input_experiments") or die "ERROR: Could not open file $input_experiments\n";
<EXPERIMENTS>;  # skip header
while (<EXPERIMENTS>) {
	chomp;
	my ($study, $policy, $experiment_id, $title, $date, $dataset_type, $library_name, $library_strategy, $library_source, $library_selection, $read_length, $instrument_model, $description) = split("\t", $_);
	
	$experiments{$experiment_id}{ID} = $experiment_id;
	$experiments{$experiment_id}{STUDY} = $study;
	$experiments{$experiment_id}{POLICY} = $policy;
	$experiments{$experiment_id}{TITLE} = $title;
	$experiments{$experiment_id}{DATE} = $date;
	$experiments{$experiment_id}{DATASET_TYPE} = $dataset_type;
	$experiments{$experiment_id}{LIBRARY_NAME} = $library_name;
	$experiments{$experiment_id}{LIBRARY_STRATEGY} = $library_strategy;
	$experiments{$experiment_id}{LIBRARY_SOURCE} = $library_source;
	$experiments{$experiment_id}{LIBRARY_SELECTION} = $library_selection;
	$experiments{$experiment_id}{READ_LENGTH} = $read_length;
	$experiments{$experiment_id}{INSTRUMENT_MODEL} = $instrument_model;
	$experiments{$experiment_id}{DESCRIPTION} = $description;
}
close(EXPERIMENTS);

# read analyses

my %analyses;
open(ANALYSES, "$input_analyses") or die "ERROR: Could not open file $input_analyses\n";
<ANALYSES>;  # skip header
while (<ANALYSES>) {
	chomp;
	my ($sample_alias, $experiment_id, $filename, $bam_md5_encrypted, $bam_md5_unencrypted, $bai_md5_encrypted, $bai_md5_unencrypted) = split("\t", $_);

	die "ERROR: Sample alias '$sample_alias' not found in file $input_samples\n" if (!exists $samples{$sample_alias});
	die "ERROR: Experiment ID '$experiment_id' not found in file $input_experiments\n" if (!exists $experiments{$experiment_id});
	
	my $analysis_id = $sample_alias."_".$experiment_id;
	$analyses{$analysis_id}{ID} = $analysis_id;
	$analyses{$analysis_id}{FILENAME} = $filename;
	$analyses{$analysis_id}{BAM_MD5_ENCRYPTED} = $bam_md5_encrypted;
	$analyses{$analysis_id}{BAM_MD5_UNENCRYPTED} = $bam_md5_unencrypted;
	$analyses{$analysis_id}{BAI_MD5_ENCRYPTED} = $bai_md5_encrypted;
	$analyses{$analysis_id}{BAI_MD5_UNENCRYPTED} = $bai_md5_unencrypted;
	$analyses{$analysis_id}{SAMPLE} = $samples{$sample_alias};
	$analyses{$analysis_id}{EXPERIMENT} = $experiments{$experiment_id};
	
	$samples{$sample_alias}{ANALYSES} = {} if (!exists $samples{$sample_alias}{ANALYSES});
	$samples{$sample_alias}{ANALYSES}{$analysis_id} = $analyses{$analysis_id};

	$experiments{$experiment_id}{ANALYSES} = {} if (!exists $experiments{$experiment_id}{ANALYSES});
	$experiments{$experiment_id}{ANALYSES}{$analysis_id} = $analyses{$analysis_id};
}
close(ANALYSES);

# read submission result file
# NOTE: this file is obtained by first submitting sample.xml and analysis.xml to EBI REST service

my $subres = "";
open(SUBRES, "$input_submission_result") or die "ERROR: Could not open file $input_submission_result\n";
while(<SUBRES>) {
	$subres .= $_;
}
close(SUBRES);

while($subres =~ /Samples: \nAccession:([^\n]+)\nBioSample Id:([^\n]+)\nAlias: ([^\n]+)\nStatus: ([^\n]+)\n/g) {
	$samples{$3}{ACCESSION} = $1;
	$samples{$3}{BIOSAMPLEID} = $2;
	$samples{$3}{STATUS} = $4;
}
while($subres =~ /Analysis: \nAccession:([^\n]+)\nAlias: ([^\n]+)\nStatus: ([^\n]+)\n/g) {
	$analyses{$2}{ACCESSION} = $1;
	$analyses{$2}{STATUS} = $3;
}


#----------------
# sample XML
#----------------

# build hash

my $schema = XML::Compile::Schema->new("/mnt/projects/hdall/results/ega-submission/schemas/SRA.sample.xsd");
$schema->importDefinitions('/mnt/projects/hdall/results/ega-submission/schemas/SRA.common.xsd');
#$schema->printIndex();
#warn $schema->template('PERL', 'SAMPLE_SET');

my $sampleSet = { seq_SAMPLE => [] };
foreach my $sample_alias (keys(%samples)) {	

	my %sample;
	$sample{SAMPLE}{alias} = $sample_alias;
	$sample{SAMPLE}{center_name} = "CCRI";
	$sample{SAMPLE}{TITLE} = $samples{$sample_alias}{TITLE};
	$sample{SAMPLE}{DESCRIPTION} = $samples{$sample_alias}{DESCRIPTION};
	
	$sample{SAMPLE}{SAMPLE_NAME} = { 
		TAXON_ID  => $samples{$sample_alias}{SUBJECT}{TAXON_ID}, 
		SCIENTIFIC_NAME => $samples{$sample_alias}{SUBJECT}{SCIENTIFIC_NAME}, 
		COMMON_NAME => $samples{$sample_alias}{SUBJECT}{COMMON_NAME}
	};
	 
	my $sampleAttributes = { seq_SAMPLE_ATTRIBUTE => [] };
	push($sampleAttributes->{seq_SAMPLE_ATTRIBUTE}, { SAMPLE_ATTRIBUTE => { TAG => "subject_id", VALUE => $samples{$sample_alias}{SUBJECT}{ID} }}); 
	push($sampleAttributes->{seq_SAMPLE_ATTRIBUTE}, { SAMPLE_ATTRIBUTE => { TAG => "phenotype", VALUE => $samples{$sample_alias}{PHENOTYPE} }}); 
	push($sampleAttributes->{seq_SAMPLE_ATTRIBUTE}, { SAMPLE_ATTRIBUTE => { TAG => "gender", VALUE => $samples{$sample_alias}{SUBJECT}{GENDER} }}); 
	push($sampleAttributes->{seq_SAMPLE_ATTRIBUTE}, { SAMPLE_ATTRIBUTE => { TAG => "cohort", VALUE => $samples{$sample_alias}{SUBJECT}{COHORT} }}); 
	push($sampleAttributes->{seq_SAMPLE_ATTRIBUTE}, { SAMPLE_ATTRIBUTE => { TAG => "age_at_diagnosis_years", VALUE => $samples{$sample_alias}{SUBJECT}{AGE_AT_DIAGNOSIS_YEARS} }}); 
	$sample{SAMPLE}{SAMPLE_ATTRIBUTES} = $sampleAttributes;
		
	push($sampleSet->{seq_SAMPLE}, \%sample); 
}
close(SAMPLES);

#print Dumper($sampleSet);

# write XML

{
	my $doc    = XML::LibXML::Document->new('1.0', 'UTF-8');
	my $write  = $schema->compile(WRITER => 'SAMPLE_SET');
	my $xml    = $write->($doc, $sampleSet);
	$xml->addChild ($doc->createAttribute ( 'xmlns:xsi' => 'http://www.w3.org/2001/XMLSchema-instance' ) );
	$xml->addChild ($doc->createAttribute ( 'xsi:noNamespaceSchemaLocation' => $schema_sample ) );
	$doc->setDocumentElement($xml);
	
	open(XML,">$output_sample_xml") or die "ERROR: Could not write to $output_sample_xml\n";
	print XML $doc->toString(1);
	close(XML);
	print STDERR "Sample XML written to $output_sample_xml\n";	
}

#----------------
# experiment XML
#----------------

# build hash

$schema = XML::Compile::Schema->new("/mnt/projects/hdall/results/ega-submission/schemas/SRA.experiment.xsd");
$schema->importDefinitions('/mnt/projects/hdall/results/ega-submission/schemas/SRA.common.xsd');
#$schema->printIndex();
#warn $schema->template('PERL', 'EXPERIMENT_SET');

my $experimentSet = { EXPERIMENT => [] };
foreach my $analysis_id (keys(%analyses)) {	
	
	my %experiment;
	$experiment{alias} = $analysis_id;
	$experiment{center_name} = "CCRI";
	$experiment{broker_name} = "EGA";
	$experiment{accession} = "";
	$experiment{TITLE} = $analyses{$analysis_id}{EXPERIMENT}{TITLE};
	$experiment{STUDY_REF} = { accession => $analyses{$analysis_id}{EXPERIMENT}{STUDY} };
	$experiment{DESIGN} = { 
		DESIGN_DESCRIPTION => $analyses{$analysis_id}{EXPERIMENT}{DESCRIPTION},
		SAMPLE_DESCRIPTOR => { refname => $analyses{$analysis_id}{SAMPLE}{ID} },
		LIBRARY_DESCRIPTOR =>  {
			LIBRARY_STRATEGY => $analyses{$analysis_id}{EXPERIMENT}{LIBRARY_STRATEGY},
			LIBRARY_SOURCE => $analyses{$analysis_id}{EXPERIMENT}{LIBRARY_SOURCE},
			LIBRARY_SELECTION => $analyses{$analysis_id}{EXPERIMENT}{LIBRARY_SELECTION},
			LIBRARY_LAYOUT => { PAIRED => {  }}
		}
	};
	$experiment{PLATFORM} = {
		ILLUMINA => { INSTRUMENT_MODEL => $analyses{$analysis_id}{EXPERIMENT}{INSTRUMENT_MODEL} }
	}; 
	$experiment{PROCESSING} = { };
	
	push($experimentSet->{EXPERIMENT}, \%experiment); 		
}

#print Dumper($experimentSet);

# write XML
{
	my $doc    = XML::LibXML::Document->new('1.0', 'UTF-8');
	my $write  = $schema->compile(WRITER => 'EXPERIMENT_SET');
	my $xml    = $write->($doc, $experimentSet);
	#$xml->addChild ($doc->createAttribute ( 'xmlns:xsi' => 'http://www.w3.org/2001/XMLSchema-instance' ) );
	$xml->addChild ($doc->createAttribute ( 'xsi:noNamespaceSchemaLocation' => $schema_experiment ) );
	$doc->setDocumentElement($xml);
	
	open(XML,">$output_experiment_xml") or die "ERROR: Could not write to $output_experiment_xml\n";
	print XML $doc->toString(1);
	close(XML);
	print STDERR "Experiment XML written to $output_experiment_xml\n";	
}

#----------------
# run XML
#----------------

# build hash

$schema = XML::Compile::Schema->new("/mnt/projects/hdall/results/ega-submission/schemas/SRA.run.xsd");
$schema->importDefinitions('/mnt/projects/hdall/results/ega-submission/schemas/SRA.common.xsd');
#$schema->printIndex();
#print $schema->template('PERL', 'RUN_SET');

my $runSet = { RUN => [] };
foreach my $analysis_id (keys(%analyses)) {	
	
	my %run;
	$run{alias} = $analysis_id;
	$run{center_name} = "CCRI";
	$run{run_center} = "CeMM";
	$run{run_date} = $analyses{$analysis_id}{EXPERIMENT}{DATE};
	$run{EXPERIMENT_REF} = { refname => $analysis_id };
	$run{RUN_TYPE}{REFERENCE_ALIGNMENT}{ASSEMBLY}{STANDARD} = { refname => "hg19" };
	$run{DATA_BLOCK}{FILES}{FILE} = { 
		filename => $analyses{$analysis_id}{FILENAME}, 
		filetype => 'bam',, 
		checksum_method => "MD5", 
		checksum => $analyses{$analysis_id}{BAM_MD5_ENCRYPTED}, 
		unencrypted_checksum => $analyses{$analysis_id}{BAM_MD5_UNENCRYPTED}
	}; 
	
	push($runSet->{RUN}, \%run); 		
}

#print Dumper($runSet);

# write XML
{
	my $doc    = XML::LibXML::Document->new('1.0', 'UTF-8');
	my $write  = $schema->compile(WRITER => 'RUN_SET');
	my $xml    = $write->($doc, $runSet);
	#$xml->addChild ($doc->createAttribute ( 'xmlns:xsi' => 'http://www.w3.org/2001/XMLSchema-instance' ) );
	$xml->addChild ($doc->createAttribute ( 'xsi:noNamespaceSchemaLocation' => $schema_run ) );
	$doc->setDocumentElement($xml);
	
	open(XML,">$output_run_xml") or die "ERROR: Could not write to $output_run_xml\n";
	print XML $doc->toString(1);
	close(XML);
	print STDERR "Run XML written to $output_run_xml\n";	
}

#----------------
# analysis XML
#----------------

# build hash

$schema = XML::Compile::Schema->new("/mnt/projects/hdall/results/ega-submission/schemas/SRA.analysis.xsd");
$schema->importDefinitions('/mnt/projects/hdall/results/ega-submission/schemas/SRA.common.xsd');
#$schema->printIndex();
#print $schema->template('PERL', 'ANALYSIS_SET');

my $analysisSet = { seq_ANALYSIS => [] };
foreach my $analysis_id (keys(%analyses)) {	
	
	my %analysis;
	$analysis{ANALYSIS}{alias} = $analysis_id;
	$analysis{ANALYSIS}{center_name} = "CCRI";
	$analysis{ANALYSIS}{broker_name} = "EGA";
	$analysis{ANALYSIS}{analysis_center} = "CeMM";
	$analysis{ANALYSIS}{analysis_date} = $analyses{$analysis_id}{EXPERIMENT}{DATE};
	$analysis{ANALYSIS}{TITLE} = $analyses{$analysis_id}{SAMPLE}{ID}." (".$analyses{$analysis_id}{EXPERIMENT}{ID}.")";
	$analysis{ANALYSIS}{DESCRIPTION} = $analyses{$analysis_id}{EXPERIMENT}{TITLE}.". ".$analyses{$analysis_id}{EXPERIMENT}{DESCRIPTION};
	$analysis{ANALYSIS}{STUDY_REF} = { accession => $analyses{$analysis_id}{EXPERIMENT}{STUDY}, refcenter => "CCRI" };
	$analysis{ANALYSIS}{SAMPLE_REF} = { refname => $analyses{$analysis_id}{SAMPLE}{ID}, refcenter => "CCRI" };
#	$analysis{ANALYSIS}{EXPERIMENT_REF} = { refname => $analysis_id };
#	$analysis{ANALYSIS}{RUN_REF} = { refname => $analysis_id };
	$analysis{ANALYSIS}{ANALYSIS_TYPE}{REFERENCE_ALIGNMENT}{ASSEMBLY}{STANDARD} = { refname => "hg19" };
	
	my @files;
	push(@files, {  
			filename => $analyses{$analysis_id}{FILENAME}, 
			filetype => 'bam', 
			checksum_method => "MD5", 
			checksum => $analyses{$analysis_id}{BAM_MD5_ENCRYPTED},
			unencrypted_checksum => $analyses{$analysis_id}{BAM_MD5_UNENCRYPTED}
	});
	if ($analyses{$analysis_id}{BAI_MD5_ENCRYPTED} ne "") {
		my $indexfile = $analyses{$analysis_id}{FILENAME};
		$indexfile =~ s/\.bam$/\.bai/;
		push(@files, { 
				filename => $indexfile, 
				filetype => 'bai', 
				checksum_method => "MD5", 
				checksum => $analyses{$analysis_id}{BAI_MD5_ENCRYPTED},
				unencrypted_checksum => $analyses{$analysis_id}{BAI_MD5_UNENCRYPTED}
		});		
	}		
	$analysis{ANALYSIS}{FILES}{FILE} = \@files;
	
	push($analysisSet->{seq_ANALYSIS}, \%analysis); 		
}

#print Dumper($analysisSet);

# write XML
{
	my $doc    = XML::LibXML::Document->new('1.0', 'UTF-8');
	my $write  = $schema->compile(WRITER => 'ANALYSIS_SET');
	my $xml    = $write->($doc, $analysisSet);
	$xml->addChild ($doc->createAttribute ( 'xmlns:xsi' => 'http://www.w3.org/2001/XMLSchema-instance' ) );
	$xml->addChild ($doc->createAttribute ( 'xsi:noNamespaceSchemaLocation' => $schema_analysis ) );
	$doc->setDocumentElement($xml);
	
	open(XML,">$output_analysis_xml") or die "ERROR: Could not write to $output_analysis_xml\n";
	print XML $doc->toString(1);
	close(XML);
	print STDERR "Analysis XML written to $output_analysis_xml\n";	
}

#----------------
# dataset XML
#----------------

$schema = XML::Compile::Schema->new("/mnt/projects/hdall/results/ega-submission/schemas/EGA.dataset.xsd");
$schema->importDefinitions('/mnt/projects/hdall/results/ega-submission/schemas/SRA.common.xsd');
$schema->printIndex();
print $schema->template('PERL', 'DATASETS');

my $datasets = { seq_DATASET => [] };
foreach my $experiment_id (keys(%experiments)) {
	
	my %dataset;
	$dataset{DATASET}{alias} = $experiment_id;
	$dataset{DATASET}{center_name} = "CCRI";
	$dataset{DATASET}{broker_name} = "EGA";
	$dataset{DATASET}{TITLE} = $experiments{$experiment_id}{TITLE};
	$dataset{DATASET}{DESCRIPTION} = $experiments{$experiment_id}{DESCRIPTION};
	$dataset{DATASET}{DATASET_TYPE} = [ $experiments{$experiment_id}{DATASET_TYPE} ];
	$dataset{DATASET}{POLICY_REF} = { accession => $experiments{$experiment_id}{POLICY}, refcenter => "CCRI" };

	my @analysis_ref;	
	foreach my $analysis_id (keys($experiments{$experiment_id}{ANALYSES})) {
		push(@analysis_ref, { accession => $analyses{$analysis_id}{ACCESSION}, refcenter => "CCRI" });
	}
	$dataset{DATASET}{ANALYSIS_REF} = \@analysis_ref;

	
	push($datasets->{seq_DATASET}, \%dataset);
}	

# write XML
{
	my $doc    = XML::LibXML::Document->new('1.0', 'UTF-8');
	my $write  = $schema->compile(WRITER => 'DATASETS');
	my $xml    = $write->($doc, $datasets);
	$xml->addChild ($doc->createAttribute ( 'xmlns:xsi' => 'http://www.w3.org/2001/XMLSchema-instance' ) );
	$xml->addChild ($doc->createAttribute ( 'xsi:noNamespaceSchemaLocation' => $schema_dataset ) );
	$doc->setDocumentElement($xml);
	
	open(XML,">$output_dataset_xml") or die "ERROR: Could not write to $output_dataset_xml\n";
	print XML $doc->toString(1);
	close(XML);
	print STDERR "Dataset XML written to $output_dataset_xml\n";	
}

#----------------
# submission XML
#----------------

# build hash

$schema = XML::Compile::Schema->new("/mnt/projects/hdall/results/ega-submission/schemas/SRA.submission.xsd");
$schema->importDefinitions('/mnt/projects/hdall/results/ega-submission/schemas/SRA.common.xsd');
#$schema->printIndex();
#warn $schema->template('PERL', 'SUBMISSION_SET');

my $submissionSet = { 
	seq_SUBMISSION => [
		{
			SUBMISSION => {
				alias => "hdall",
				center_name => "CCRI",
				broker_name => "EGA",
				ACTIONS => {
					seq_ACTION => [
#						{
#							ACTION => {
#								ADD => {
#									source => "sample.xml", schema => "sample"
#								}
#							}
#						},

#						{
#							ACTION => {
#								VALIDATE => {
#									source => "experiment.xml", schema => "experiment"
#								}
#							}
#						},

#						{
#							ACTION => {
#								VALIDATE => {
#									source => "run.xml", schema => "run"
#								}
#							}
#						},

#						{
#							ACTION => {
#								ADD => {
#									source => "analysis.xml", schema => "analysis"
#								}
#							}
#						},

						{
							ACTION => {
								ADD => {
									source => "dataset.xml", schema => "dataset"
								}
							}
						},
						{
							ACTION => {
								PROTECT => { }
							}
						}
					]
				}
			}
		}
	] 
};

# write XML

{
	my $doc    = XML::LibXML::Document->new('1.0', 'UTF-8');
	my $write  = $schema->compile(WRITER => 'SUBMISSION_SET');
	my $xml    = $write->($doc, $submissionSet);
	$xml->addChild ($doc->createAttribute ( 'xmlns:xsi' => 'http://www.w3.org/2001/XMLSchema-instance' ) );
	$xml->addChild ($doc->createAttribute ( 'xsi:noNamespaceSchemaLocation' => $schema_submission ) );
	$doc->setDocumentElement($xml);
	
	open(XML,">$output_submission_xml") or die "ERROR: Could not write to $output_submission_xml\n";
	print XML $doc->toString(1);
	close(XML);
	print STDERR "Submission XML written to $output_submission_xml\n";
}
