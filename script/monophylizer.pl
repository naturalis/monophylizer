#!/usr/bin/perl
use strict;
use warnings;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use Pod::Usage;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# release version for ExtUtils::MakeMaker to parse
our $VERSION = '0.1';

=head1 NAME

monophylizer.pl - assesses taxonomic monophyly on Barcode of Life trees

=head1 SYNOPSIS

 monophylizer.pl [options]
 
=head1 OPTIONS AND ARGUMENTS

=over

=item B<-infile> <file>

A tree file, usually in Newick format. Required.

=item B<-format> <newick|nexus|nexml|phyloxml>

Optional argument to specify the tree file format. By default the Newick format is used.

=item B<-metadata> <file>

Optional argument to provide the location of a tab-separated spreadsheet with per-taxon 
metadata.

=item B<-separator> <character>

Optional argument to specific the character that separates the taxon name from any
additional metdata (such as sequence IDs) in leaf labels. By default this is the
pipe symbol: '|'.

=item B<-comments>

Optional flag to treat square brackets as opaque strings, not comments, default: true.

=item B<-quotes>

Optional flag to treat quotes as opaque strings, default: true.

=item B<-whitespace> 

Optional flag to treat whitespace as opaque strings, default: true.

=item B<-trinomials>

Optional flag to include subspecific epithets in taxa, default: false.

=item B<-astsv>

Optional flag to set output as TSV regardless whether running as CGI, default: false.
This is only available when running under CGI.

=item B<-verbose>

Influences how verbose the script is. By default, only warning messages are emitted.
When this flag is used once, also informational messages are emitted. When used twice,
also debugging messages.

=item B<-help>

Prints usage message and quits, only available when running on command line.

=back

=head1 DESCRIPTION

This script assesses whether the species in a phylogenetic tree are monophyletic. It 
reports the output in a table that lists for each species its status (monophyletic,
polyphyletic or paraphyletic), and if not monophyletic, which other species it is
entangled with. In addition, sequence identifiers and any other input metadata are
reported. How the assessment is made algorithmically is described below, in the section
on the C<do_assessment> subroutine.

=head1 SUBROUTINES

This section describes the various subroutines that implement the functionality of the
script. This information is generally only of interest to developers.

=over

=item main

The main subroutine takes the following steps:

 - process_args:     get the input from CGI or command line
 - make_taxa:        create a taxa block from the input tree's leaf labels
 - read_spreadsheet: join the taxa with additional metadata, if any
 - index_nodes:      applies pre- and post-order labels, aggregates "stop nodes"
 - do_assessment:    assess whether stop node distribution implies mono/para/poly
 - [print_html:      print output as HTML table (only under CGI without conneg)]
 - [or print a tab-separated spreadsheet]

=cut

sub main {

	# process the arguments, either from the command line
	# or from CGI parameters
	my %args = process_args();

	# create and populate a Bio::Phylo::Taxa object by
	# parsing tip labels in the provided tree
	my $taxa = make_taxa(@args{qw(tree trinomials separator factory)});

	# read additional metadata, if any
	read_spreadsheet($taxa,$args{'metafh'});

	# do the tree traversal to identify "stopnodes"
	index_nodes($args{'tree'});

	# do the final assessment based on the topology of
	# the stopnodes. assemble the list of tanglees
	# accordingly.
	my @result = do_assessment($taxa,$args{'tree'});

	# HTML is printed when running as a web service
	if ( $args{'cgi'} and not $args{'astsv'} ) {
		print_html(@result);
	}
	else {

		# this header is printed when we are running as a web service
		# and the client has requested TSV output, either by clicking
		# a checkbox (i.e. passing 'astsv=checked' in the query string)
		# or by explicitly asking for 'text/plain' in the Accept header
		# (i.e. through content negotiation).
		if ( $args{'cgi'} ) {
			print "Content-type: text/plain\n\n";
		}
	
		# the remainder is the same regardless whether we are in the
		# shell or on a server
		my @header = qw(Species Assessment Tanglees IDs Metadata);
		print join("\t",@header), "\n";	
		for my $r ( sort { $a->[0] cmp $b->[0] } @result ) {
			print join("\t",@$r), "\n";
		}
	}
}
main();

=item process_args

Processes command line or CGI arguments, returns a command object with argument fields. 
The script determines whether or not it is running under CGI by checking whether the 
hidden form field 'cgi' is set to a true value. The HTML page that composes the CGI 
request turns this flag on, so any other clients (e.g. cURL on the command line) need to 
ensure they do the same.

[Optional] arguments read from command line or CGI:

 infile:      a tree file
 [metadata:   a tab-separated spreadsheet with per-taxon metadata]
 [format:     tree file format, default: 'newick']
 [verbose:    verbosity level, default: WARN]
 [separator:  symbol that separates leaf labels from other fields, default: '|']
 [comments:   flag to treat square brackets as opaque strings, not comments, default: true]
 [quotes:     flag to treat quotes as opaque strings, default: true]
 [whitespace: flag to treat whitespace as opaque strings, default: true]
 [trinomials: flag to include subspecific epithets in taxa, default: false]
 [astsv:      emit output as TSV regardless whether running as CGI, default: false]
 [help:       prints usage message and quits, only available when running on command line]

Returned command object fields:

 tree:       a Bio::Phylo::Forest::Tree object
 log:        a Bio::Phylo::Util::Logger object
 factory:    a Bio::Phylo::Factory object 
 metafh:     a file handle for the metadata spreadsheet, if any
 separator:  symbol that separates leaf labels from other fields
 trinomials: flag to include subspecific epithets in taxa
 cgi:        flag to indicate whether we are running as CGI
 astsv:      emit output as TSV regardless whether running as CGI

=cut

sub process_args {
	# this will be set to true if we are running in a CGI environment
	# so that we then correctly set the MIME-type in the output
	my $as_cgi_process;

	# process command line arguments
	my $verbosity = WARN;
	my $format = 'newick';
	my $separator = '|';

	# treat comment characters as opaque: BOLD does not
	# introduce [newick comments], but it *does*
	# create records that potentially have square
	# brackets in them, which need to be treated as
	# any other character without hidden semantics
	my $comments = 1;

	# likewise: treat quotes as opaque
	my $quotes = 1;

	# treat whitespace as part of the name or metadata:
	# BOLD does not introduce pretty indentation in the 
	# newick string (e.g. as in the NCBI classification 
	# tree) *and* it does not escape whitespace by quoting,
	# so we need to retain it
	my $whitespace = 1;

	# optionally, the user may supply additional metadata
	# in tab-separated format. this needs to be cross-referenced
	# with species names in the final result
	my $metadata;

	# this is the input tree file, by default this will be
	# a file with a single newick tree in it that uses the
	# formatting variant implemented by BOLD
	my $infile;

	# this flag indicates whether to treat trinomials
	# (i.e. species + subspecific epithet) as separate
	# taxa. by default we do not do this: anything after
	# putative genus + species is ignore
	my $trinomials = 0;
	
	# this flag indicates whether to write out tab-separated
	# data in web service requests. This is handy for clients
	# wanting to import the results directly into a spreadsheet
	# program. This flag is either toggled by an HTML checkbox
	# (i.e. a query parameter) or using content negotiations
	my $astsv;

	my $help;
	GetOptions(
		'infile=s'    => \$infile,
		'verbose+'    => \$verbosity,
		'format=s'    => \$format,
		'separator=s' => \$separator,
		'comments'    => \$comments,
		'whitespace'  => \$whitespace,
		'trinomials'  => \$trinomials,
		'astsv'       => \$astsv,
		'metadata=s'  => \$metadata,
		'help|?'      => \$help,
		'quotes'      => \$quotes,
	);

	# input file handle for tree and metadata
	my $infh;
	my $metafh;

	# process CGI arguments
	my $cgi = CGI->new;
	if ( $as_cgi_process = $cgi->param('cgi') ) {
		
		# this will spit out results in error_log. 
		# let's keep this to a minimum
		$verbosity = WARN;
		
		# first check the header and assign that way
		if ( CGI->can('Accept') ) {
			$astsv = ( CGI::Accept('text/plain') == 1.0 );
		}
		elsif ( CGI->can('CGI::accept') ) {
			$astsv = ( CGI::accept('text/plain') == 1.0 );			
		}	
			
		# retrieve simple parameter values
		$format     = $cgi->param('format');
		$separator  = $cgi->param('separator');
		$comments   = $cgi->param('comments');
		$whitespace = $cgi->param('whitespace');
		$trinomials = $cgi->param('trinomials');
		$infile     = $cgi->param('infile');
		$metadata   = $cgi->param('metadata');
		$quotes     = $cgi->param('quotes');
		
		# if people compose their own request, completely omitting the
		# parameter, we would still have the conneg one because this
		# one would be undefined. if they uncheck the box it would be
		# the empty string.
		$astsv = $cgi->param('astsv') if defined $cgi->param('astsv');
		
		# TO DO: maybe make it so that users can also provide a location
		# for the input data, so that the monophylizer can be embedded
		# in BioVeL/Taverna workflows.
		$infh       = $cgi->upload('infile')->handle;
		$metafh     = $cgi->upload('metadata')->handle if $metadata;
	}
	else {
		open $infh,   '<', $infile   or die $!;
		open $metafh, '<', $metadata or die $! if $metadata;
	}

	# emit help message if run with --help, -help, -h or -?
	if ( $help ) {
		print <<"USAGE";
	$0 -i <tree file> [-f <format>] [-s <separator>] [-m <metadata>] [--verbose]

	Default file format is newick, default separator between species name and identifier
	is the '|' symbol.

USAGE
		exit 0;
	}

	# instantiate helper objects
	my $fac = Bio::Phylo::Factory->new;
	my $log = Bio::Phylo::Util::Logger->new(
		'-level' => $verbosity,
		'-class' => [ 'main' ],
	);
	$log->info("going to read $format tree from $infile");
	my $tree = parse_tree(
		'-format'          => $format,
		'-handle'          => $infh,
		'-as_project'      => 1,
		'-ignore_comments' => !!$comments,
		'-ignore_quotes'   => !!$quotes,
		'-keep_whitespace' => !!$whitespace,
	);
	$log->info("done reading tree: $tree");
	return 
		'tree'       => $tree,
		'log'        => $log,
		'metafh'     => $metafh,
		'separator'  => $separator,
		'trinomials' => $trinomials,
		'factory'    => $fac,
		'cgi'        => $as_cgi_process,
		'astsv'      => $astsv;
}

=item make_taxa

Given a L<Bio::Phylo::Forest::Tree> object, makes a L<Bio::Phylo::Taxa> object that 
contains all the distinct species names, with their ID annotations. The species names
are extracted from the leaf labels by splitting these on a $separator token and taking
the first segment, then splitting that on underscores or spaces to get the taxonomic
name parts (e.g. genus, species, subspecies). If the $trinomials flag is set to true,
three parts of the name are concatenated, otherwise two. For each distinct name that
is created like this, a corresponding L<Bio::Phylo::Taxa::Taxon> object is created, and
each node whose label contains this name is linked to it. The taxon object aggregates
all the sequence IDs (i.e. the part after $separator) into an array that it stores in
the 'ids' field.

Arguments:

 $tree:       input tree
 $trinomials: if true, subspecific epithets are included in the taxon name
 $separator:  the token on which to split the leaf label string
 $fac:        factory object that creates taxa and taxon objects

=cut

sub make_taxa {
	my ( $tree, $trinomials, $separator, $fac ) = @_;
	my %taxa;
	for my $tip ( @{ $tree->get_terminals } ) {
		my $label = $tip->get_name;
		if ( $label =~ /^'?(.+?)\Q$separator\E(.+?)'?$/ ) {
			my ( $name, $meta ) = ( $1, $2 );
			
			# sometimes it's underscores, sometimes it's spaces
			my @parts = split /(?:_|\s)/, $name;
			
			# make a bi- or trinomial taxon name
			my $taxon_name;
			if ( $trinomials and scalar @parts > 2 ) {
				$taxon_name = join ' ', @parts[0..2];
			}
			else {
				no warnings 'uninitialized';
				$taxon_name = join ' ', @parts[0,1];
				if ( not $parts[0] or not $parts[1] ) {
					warn "not a binomial: @parts";
				}
			}
			
			# fetch or instantiate the taxon by that name
			my $taxon = $taxa{$taxon_name} //= $fac->create_taxon('-name'=>$taxon_name);
			
			# attach the sequence ID
			my $ids = $taxon->get_generic('ids') || [];
			push @$ids, $meta;
			$taxon->set_generic( 'ids' => $ids );
			
			# link node to taxon
			$tip->set_taxon($taxon);
		}
	}
	my $taxa = $fac->create_taxa;
	$taxa->insert($_) for values %taxa;
	return $taxa;
}

=item read_spreadsheet

Given a L<Bio::Phylo::Taxa> object and a file handle, reads the handle as a tab-separated
spreadsheet (no header) whose first field is a name in the taxa object. Attaches all 
subsequent fields to the corresponding taxon.

Arguments:
 
 $taxa:   a Bio::Phylo::Taxa object
 $metafh: a file handle of a tab-separated spreadsheet

=cut

sub read_spreadsheet {
	my ( $taxa, $metafh ) = @_;
	if ( $metafh ) {
		my $line = 1;
		while(<$metafh>) {
			chomp;
			my @record = split /\t/, $_;
			my $name = shift @record;
			if ( my $taxon = $taxa->get_by_name($name) ) {
				$taxon->set_generic( 'meta' => \@record );
			}
			else {
				die "Can't locate '$name' on line $line of the metadata in the tree"; 
			}
			$line++;
		}
	}	
}

=item index_nodes

Given a tree, applies pre- and post-order indices and aggregates for each taxon all the
nodes where it coalesces with at least one other species (so-called "stop nodes", cf.
L<http://biophylo.blogspot.nl/2013/04/algorithm-for-distinguishing-polyphyly.html>).

=cut

sub index_nodes {
	my $tree = shift;
	my $counter = 1;
	$tree->visit_depth_first(
		'-pre' => sub { 
			
			# apply left index
			shift->set_generic( 'left'  => $counter++ ) 
		},
		'-pre_sister' => sub {
		
			# apply right index
			shift->set_generic( 'right' => $counter++ );
		},
		'-post' => sub { 
			my $node = shift;
			$node->set_generic( 'right' => $counter++ ) if not $node->get_generic('right');
			
			# the node is terminal
			if ( my $taxon = $node->get_taxon ) {
			
				# store its taxon by index
				$node->set_generic( 'taxa' => { $taxon->get_id => $taxon } );
			}
			else {
			
				# merge all taxa from the immediate children
				my %merged;
				for my $c ( @{ $node->get_children } ) {
					my %taxa = %{ $c->get_generic('taxa') };
					$merged{$_} = $taxa{$_} for keys %taxa;
				}
				$node->set_generic( 'taxa' => \%merged );
								
				# if two or more distinct taxa, it's a 'stopnode':
				# http://biophylo.blogspot.nl/2013/04/algorithm-for-distinguishing-polyphyly.html
				if ( scalar(keys(%merged)) > 1 ) {

					# tell all the subtended taxa about their stopnode				
					for my $taxon ( values %merged ) {
						my $stopnodes = $taxon->get_generic('stopnodes') || [];
						push @$stopnodes, $node;
						$taxon->set_generic( 'stopnodes' => $stopnodes );
					}
				}				
			}
		},
	);
}

=item do_assessment

Given a L<Bio::Phylo::Taxa> object and a L<Bio::Phylo::Forest::Tree> object, iterates
over all taxa. For each taxon, all the leaf nodes that belong to that taxon are 
collected. For this set of leaf nodes, assesses whether these form a monophyletic 
group. This is done by finding the MRCA of these leaf nodes, then fetch all descendants
of the MRCA. If the set of descendants is the same side as the set of leaf nodes, the
assessment for this taxon is C<monophyletic>.

If the taxon is not monophyletic, the set of "stop nodes" is retrieved. These are all
the internal nodes where the focal taxon coalesces with at least one other taxon, sorted
in a post-order traversal. It then assesses for each stop node whether it descends from
the next stop node in the sorted list. This is done by checking that the left (pre-order)
index of the focal node is larger than that of the next, and the right (post-order) index
of the focal node is smaller than the next. During this testing iteration, the stop nodes
are binned in distinct paths from the tips to the root. If there's more than one distinct 
bin/path, the taxon is considered polyphyletic, otherwise paraphyletic.

For non-monophyletic taxa, the final step is to then collect all the other taxa with which
the focal taxon is entangled. It does this by fetching for each bin/path the first (most
recent) node and taking its subtended taxa. The union of these sets of subtended taxa 
forms the set of 'tanglees'.

Arguments:

 $taxa: a Bio::Phylo::Taxa object
 $tree: a Bio::Phylo::Forest::Tree object

Returns a two-dimensional array (i.e. a table) where each record consists of the fields:

 - name of the focal taxon
 - status, i.e. 'monophyletic', 'paraphyletic' or 'polyphyletic'
 - comma separated string with 'tanglees', if any (otherwise empty string)
 - comma separated string with sequence IDs
 - comma separated string with additional metadata

=cut

sub do_assessment {
	my ($taxa,$tree) = @_;
	my @result;
	
	# iterate over taxa
	for my $taxon ( @{ $taxa->get_entities } ) {
		my $name  = $taxon->get_name;
		my $tips  =	$taxon->get_nodes; # returns all leaf nodes that map to focal taxon
		my $tid   = $taxon->get_id;
		
		# makes CSV strings of sequence IDs and additional metadata
		my $ids = join ',', @{ $taxon->get_generic('ids') };
		my $meta;
		if ( $taxon->get_generic('meta') ) {
			$meta = join ',', @{ $taxon->get_generic('meta') };
		}
		else {
			$meta = '';
		}
		
		# simplest case: taxon is monophyletic
		if ( scalar(@{$tips}) == 1 || scalar(@{$tips}) == scalar(@{$tree->get_mrca($tips)->get_descendants}) ) {
			push @result, [ $name, 'monophyletic', '', $ids, $meta ];
		}
		
		# taxon is NOT monophyletic
		else {
			my @sn = @{ $taxon->get_generic('stopnodes') };
			my @bins = ( [] );
			
			# here we traverse the array of stopnodes from more recent
			# to more ancestral and bin them by path
			for my $i ( 0 .. $#sn - 1 ) {
				my $l1 = $sn[$i]->get_generic('left');
				my $r1 = $sn[$i]->get_generic('right');
				my $l2 = $sn[$i+1]->get_generic('left');
				my $r2 = $sn[$i+1]->get_generic('right');
				
				# test if the focal one indeed nests within the next one
				if ( $l1 > $l2 and $r1 < $r2 ) {
					push @{ $bins[-1] }, $sn[$i];
				}
				else {
				
					# the focal node does NOT nest within the next one,
					# we need to start a new bin if the previous one 
					# already contains things
					if ( @{ $bins[-1] } ) {
						push @bins, [ $sn[$i] ];
					}
					else {
						$bins[-1]->[0] = $sn[$i];
					}
				}
			}
			
			# if there is only one bin: para, otherwise: poly
			my $status = scalar(@bins) == 1 ? 'paraphyletic' : 'polyphyletic';
			
			# now build the set of tanglees
			my %tanglees;
			for my $bin ( @bins ) {
				my %taxa = %{ $bin->[0]->get_generic('taxa') };
				$tanglees{$_} = $taxa{$_} for keys %taxa;
			}
			delete $tanglees{$tid};
			my @tanglees = map { $_->get_name } values %tanglees;
			push @result, [ $name, $status, join(',',@tanglees), $ids, $meta ];
		}
	}
	return @result;
}

=item print_html

Prints the result as an HTML page with table. The page expects that the JavaScript
library "sorttable.js", which allows sorting on table columns, is located in the
document root of the web server. The table has the following columns:

 - name of the focal taxon
 - status, i.e. 'monophyletic', 'paraphyletic' or 'polyphyletic'
 - comma separated string with 'tanglees', if any (otherwise empty string)
 - comma separated string with sequence IDs
 - comma separated string with additional metadata

=cut

sub print_html {
	my @result = @_;
	print <<"HEADER";
Content-type: text/html

<html>
	<head>
		<title>Monophyly assessment result</title>
		<style type="text/css">
			td, th { whitespace: nowrap; font-size: small }
			body { font-family: verdana, arial, sans-serif; font-size: small }
			th { text-align: left }
		</style>
		<script type="text/javascript" src="/sorttable.js"></script>
	</head>
	<body>
		<h1>Results, March 2014 algorithm version  (click column headers to sort)</h1>
		<table class="sortable"><tr>
HEADER
	for my $column ( qw(Species Assessment Tanglees IDs Metadata) ) {
		print "<th>$column</th>";
	}
	print '</tr>', "\n";
	for my $row ( sort { $a->[0] cmp $b->[0] } @result ) {
		print '<tr>', map { "<td>$_</td>" } @$row;
		print "</tr>\n";
	}
	print "</table></body></html>"	
}

=back

=cut
