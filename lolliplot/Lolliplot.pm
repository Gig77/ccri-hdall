#----------------------------------
# $Authors: dlarson bshore $
# $Date: 2008-09-16 16:33:54 -0500 (Tue, 16 Sep 2008) $
# $Revision: 38655 $
# $URL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/MG/MutationDiagram.pm $
#----------------------------------
package Lolliplot;

#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
use strict;
#use warnings;
no warnings;
use Carp;

use lib "$ENV{HOME}/git/gms-core/lib/perl";

use FileHandle;
use Genome;

use SVG;
use Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::View;
use Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Backbone;
use Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Domain;
use Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Mutation;
use Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Legend;
use Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::LayoutManager;

#------------------------------------------------
sub new {
    my ($class, %arg) = @_;

    my $self = {
        _basename => $arg{basename} || '',
        _output_suffix => $arg{suffix} || '',
        _transcripts => $arg{transcripts} || '',
        _output_directory => $arg{output_directory} || '.',
        _vep_frequency_field => $arg{vep_frequency_field},
        _max_display_freq => $arg{max_display_freq},
        _lolli_shape => $arg{lolli_shape},
        _domains => $arg{domains},
        _variants => $arg{variants},
    };

    unless ($self->{_transcripts}){
        confess "No reference transcripts supplied, abort!";
    }

    my @hugos = ();
    if (defined($arg{hugos})) {
        @hugos = split(',',$arg{hugos});
    }
    unless (scalar(@hugos)) {
        @hugos = qw( ALL );
    }
    $self->{_hugos} = \@hugos;
    bless($self, ref($class) || $class);
    
    die "No mutation file passed to $class" unless $arg{variants};

    $self->_parse_variants();
    $self->MakeDiagrams();
    return $self;
}

sub _get_domains {
    my ($self, $transcript_name) = @_;

    if (!defined $self->{_transcripts}{$transcript_name}) {
        warn "No transcript found for $transcript_name";
        return;
    }

    my @domains;
    for my $dk (keys(%{$self->{_domains}{$transcript_name}})) 
    {
        push @domains, {
            name => $self->{_domains}{$transcript_name}{$dk}{name},
            source => $self->{_domains}{$transcript_name}{$dk}{source},
            start => $self->{_domains}{$transcript_name}{$dk}{start},
            end => $self->{_domains}{$transcript_name}{$dk}{end}
        };
    }

    return @domains;
}

sub _add_mutation {
    my ($self, %params) = @_;
    $self->{_data} = {} unless defined $self->{_data};
    my $data = $self->{_data};

    my $hugo = $params{hugo};
    my $transcript_name = $params{transcript_name};
    my $protein_length = $params{protein_length};
    my $protein_position = $params{protein_position};
    my $mutation = $params{mutation};
    my $class = $params{class};
    my $domains = $params{domains};
    my $frequency = $params{frequency} || 1;

    print "Adding mutation HUGO: $hugo TRANSCRIPT: $transcript_name POS: $protein_position MUTATION: $mutation CLASS: $class PROT_LENGTH: $protein_length\n";

    $data->{$hugo}{$transcript_name}{length} = $protein_length;
    push @{$data->{$hugo}{$transcript_name}{domains}}, @$domains;

    if (defined($protein_position)) {
        unless (exists($data->{$hugo}{$transcript_name}{mutations}{$mutation})) {
            $data->{$hugo}{$transcript_name}{mutations}{$mutation} =
            {
                res_start => $protein_position,
                class => $class,
            };
        }
        $data->{$hugo}{$transcript_name}{mutations}{$mutation}{frequency} += $frequency;
    }
}

sub argmin(@) { # really perl?
    my @arr = @_;
    return unless @arr;
    return 0 if @arr <= 1;

    my $minidx = 0;
    for my $i (1..$#arr) {
        $minidx = $i if $arr[$i] < $arr[$minidx];
    }
    return $minidx;
}

sub _get_vep_mutation_class {
    my $type = shift;
	
	return 'ESSENTIAL_SPLICE_SITE' if ($type =~ /^(SPLICE_SITE_ACCEPTOR|SPLICE_SITE_DONOR)$/);
	return 'FRAMESHIFT_CODING' if ($type =~ /^(FRAME_SHIFT)$/);
	return 'STOP_GAINED' if ($type =~ /^(STOP_GAINED)$/);
	return 'IN_FRAME_INS' if ($type =~ /^(CODON_CHANGE_PLUS_CODON_INSERTION|CODON_INSERTION)$/);
	return 'IN_FRAME_DEL' if ($type =~ /^(CODON_CHANGE_PLUS_CODON_DELETION|CODON_DELETION)$/);
	return 'NON_SYNONYMOUS_CODING' if ($type =~ /^(NON_SYNONYMOUS_CODING|NON_SYNONYMOUS_START|START_LOST)$/);
	return 'SYNONYMOUS_CODING' if ($type =~ /^(SYNONYMOUS_START|SYNONYMOUS_CODING|SYNONYMOUS_STOP)$/);
	
	return undef;	
}

sub _get_vep_extra_fields_hash {
    my $extra = shift;
    return { map { split("=") } split(";", $extra) }
}

sub _parse_variants {
    my $self = shift;

    my $graph_all = $self->{_hugos}->[0] eq 'ALL' ? 1 : 0;
    my %hugos;
    unless($graph_all) {
        %hugos = map {$_ => 1} @{$self->{_hugos}}; #convert array to hashset
    }

	foreach my $transcript_name (keys(%{$self->{_transcripts}}))
	{
		my $hugo = $self->{_transcripts}{$transcript_name}{hugo};
		next if (!$hugo or (!$graph_all and !exists($hugos{$hugo})));
		
		foreach my $v (keys(%{$self->{_variants}{$transcript_name}}))
		{
			my ($patient, $var_type, $status, $chr, $pos, $aa_change) = split("\t", $v);
			my ($effect) = split("\t", $self->{_variants}{$transcript_name}{$v});
			my ($orig_aa, $protein_pos, $new_aa) = $aa_change =~ /([^\d]?)(\d+)(.?)/;			
			my $class = _get_vep_mutation_class($effect); 
			die "Could not map mutation class: $effect\n" if (!$class);
	        my @domains = $self->_get_domains($transcript_name);

	        #add to the data hash for later graphing
            $orig_aa |= '';
            $new_aa |= '';
            my $mutation = join($protein_pos, $orig_aa, $new_aa);
            next if $mutation eq '--';
            my $extra = ''; # _get_vep_extra_fields_hash($fields[-1]);
            my $frequency = 1;
            $frequency = 1; # $extra->{$self->{_vep_frequency_field}} if exists $extra->{$self->{_vep_frequency_field}};
           next if !$frequency;

            $self->_add_mutation(
               	hugo => $hugo,
               	transcript_name => $self->{_transcripts}{$transcript_name}{'refseq'} ?  $self->{_transcripts}{$transcript_name}{'refseq'} : $transcript_name,
               	protein_length => $self->{_transcripts}{$transcript_name}{protlen},
               	protein_position => $protein_pos,
               	mutation => $mutation,
               	class => $class,
               	domains => \@domains,
               	frequency => $frequency
           	);
    	}
	}
}

sub Data {
    my ($self) = @_;
    return $self->{_data};
}

sub MakeDiagrams {
    my ($self) = @_;
    my $data = $self->{_data};
    my $basename = join("/", $self->{_output_directory}, $self->{_basename});
    my $suffix = $self->{_output_suffix};
    foreach my $hugo (keys %{$data}) {
        foreach my $transcript (keys %{$data->{$hugo}}) {
            unless($self->{_data}{$hugo}{$transcript}{length}) {
                warn "$transcript has no protein length and is likely non-coding. Skipping...\n";
                next;
            }
            my $svg_file = $basename . $hugo . '_' . $transcript . "$suffix.svg";
            my $svg_fh = new FileHandle;
            unless ($svg_fh->open (">$svg_file")) {
                die "Could not create file '$svg_file' for writing $$";
            }
            $self->Draw($svg_fh,
                $hugo, $transcript,
                $self->{_data}{$hugo}{$transcript}{length},
                $self->{_data}{$hugo}{$transcript}{domains},
                $self->{_data}{$hugo}{$transcript}{mutations}
            );
            $svg_fh->close();
        }
    }
    return $self;
}

sub Draw {
    my ($self, $svg_fh, $hugo, $transcript, $length, $domains, $mutations) = @_;
    $DB::single = 1;
    my $document = Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::View->new(width=>'800',height=>'600',
        'viewport' => {x => 0, y => 0,
            width => 800,
            height => 600},
        left_margin => 50,
        right_margin => 50,
        id => "main_document");
    my $svg = $document->svg;

    my $backbone = Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Backbone->new(parent => $document,
        gene => $hugo,
        protein_length => $length,
        backbone_height => 50,
        style => {fill => 'none', stroke => 'black'},
        id => "protein_diagram",
        $document->content_view);
    $backbone->draw;

    my @colors = qw( aliceblue azure blanchedalmond burlywood coral cyan darkgray darkmagenta darkred darkslategray deeppink dodgerblue fuchsia goldenrod grey indigo lavenderblush lightcoral lightgreen lightseagreen lightsteelblue mediumblue mediumslateblue midnightblue olivedrab palegoldenrod papayawhip plum rosybrown sandybrown slategrey tan );
    my $color = 0;
    my %domains;
    my %domains_location;
    my %domain_legend;
    foreach my $domain (sort {sort_domain($a,$b)}  @{$domains}) {
        if ($domain->{source} eq 'superfamily') {
            next;
        }
        my $domain_color;
        if (exists($domain_legend{$domain->{name}})) {
            $domain_color = $domain_legend{$domain->{name}};
        } else {
            if($color == @colors) {
                #protect against an array overrun
                $color = 0;
            }
            $domain_color = $colors[$color++];
            $domain_legend{$domain->{name}} = $domain_color;
        }
        if (exists($domains_location{$domain->{name}}{$domain->{start} . $domain->{end}})) {
            next;
        }
        $domains_location{$domain->{name}}{$domain->{start} . $domain->{end}} += 1;
        $domains{$domain->{name}} += 1;
        my $subid = '';
        if ($domains{$domain->{name}} > 1) {
            $subid = '_subid' . $domains{$domain->{name}};
        }
        my $test_domain = Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Domain->new(backbone => $backbone,
            start_aa => $domain->{start},
            stop_aa => $domain->{end},
            id => 'domain_' . $domain->{name} . $subid,
            text => $domain->{name},
            style => { fill => $domain_color,
                stroke => 'black'});
        $color++;
        $test_domain->draw;
    }
    my $domain_legend =
    Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Legend->new(backbone => $backbone,
        id => 'domain_legend',
        x => $length / 2,
        values => \%domain_legend,
        object => 'rectangle',
        style => {stroke => 'black', fill => 'none'});
    $domain_legend->draw;

    my @mutation_objects;
    my %mutation_class_colors = (
        # tgi annotator colors
        'frame_shift_del' => 'darkolivegreen',
        'frame_shift_ins' => 'crimson',
        'in_frame_del' => 'gold',
        'in_frame_ins' => 'darkgray',
        'missense' => 'cornflowerblue',
        'nonsense' => 'goldenrod',
        'splice_site_del' => 'orchid',
        'splice_site_ins' => 'saddlebrown',
        'splice_site_snp' => 'lightpink',

        # vep annotator colors
        'essential_splice_site' => 'orchid',
        'frameshift_coding' => 'darkolivegreen',
        'stop_gained' => 'goldenrod',
        'non_synonymous_coding' => 'cornflowerblue',
        'synonymous_coding' => 'white',


        'other' => 'black',
    );
    my %mutation_legend;
    my $max_frequency = 0;
    my $max_freq_mut;
    foreach my $mutation (keys %{ $mutations}) {
        $mutations->{$mutation}{res_start} ||= 0;
        my $mutation_color = $mutation_class_colors{lc($mutations->{$mutation}{class})};
        $mutation_color ||= $mutation_class_colors{'other'};
        $mutation_legend{$mutations->{$mutation}{class}} = $mutation_color;
        my $mutation_element =
        Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Mutation->new(backbone => $backbone,
            id => $mutation,
            start_aa => $mutations->{$mutation}{res_start},
            text => $mutation,
            frequency => $mutations->{$mutation}{frequency},
            color => $mutation_color,
            style => {stroke => 'black', fill => 'none'},
            max_freq => $self->{_max_display_freq},
            shape => $self->{_lolli_shape},
        );


        #jitter labels as a test
        push @mutation_objects, $mutation_element;
        if($mutations->{$mutation}{frequency} > $max_frequency) {
            $max_frequency = $mutations->{$mutation}{frequency};
            $max_freq_mut = $mutation_element;
        }
    }
    map {$_->vertically_align_to($max_freq_mut)} @mutation_objects;
    my $mutation_legend =
    Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::Legend->new(backbone => $backbone,
        id => 'mutation_legend',
        x => 0,
        values => \%mutation_legend,
        object => 'circle',
        style => {stroke => 'black', fill => 'none'});
    $mutation_legend->draw;


    my $layout_manager = Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::LayoutManager->new(iterations => 1000,
        max_distance => 13, spring_constant => 6, spring_force => 1, attractive_weight => 5 );
    $layout_manager->layout(@mutation_objects);

    map {$_->draw;} (@mutation_objects);

    # now render the SVG object, implicitly use svg namespace
    print $svg_fh $svg->xmlify;
}

sub domain_length {
    my $domain = shift;
    return $domain->{end} - $domain->{start} + 1;
}

sub sort_domain {
    my ($a, $b) = @_;
    my $sort_val = domain_length($b) <=> domain_length($a);
    if($sort_val == 0) {
        return $b->{name} eq $a->{name};
    }
    else {
        return $sort_val;
    }
}
