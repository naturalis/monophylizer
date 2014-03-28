#!/usr/bin/perl
use strict;
use warnings;
use CGI;
use Data::Dumper;
use CGI::Carp 'fatalsToBrowser';
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

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

# does the poly/para assessment as per:
# http://biophylo.blogspot.nl/2013/04/algorithm-for-distinguishing-polyphyly.html
sub do_assessment {
	my ($taxa,$tree) = @_;
	my @result;
	$taxa->visit(sub{
		my $taxon  = shift;
		my $name   = $taxon->get_name;
		my $status = 'monophyletic';
		my @ids    = @{ $taxon->get_generic('ids') };
		my @meta   = @{ $taxon->get_generic('meta') } if $taxon->get_generic('meta');
		my @tanglees;	
		if ( my $nodes = $taxon->get_generic('stopnodes') ) {
			if ( scalar(@$nodes) > 1 ) {
				my @n = sort { $a->get_generic('left') <=> $b->get_generic('left') } @$nodes;
				$status = 'paraphyletic';
				for my $i ( 0 .. $#n - 1 ) {
					$status = 'polyphyletic' if $n[$i]->get_generic('right') < $n[$i+1]->get_generic('right');
				}
				my $mrca = $tree->get_mrca( $taxon->get_nodes );
				my %t = map { $_->get_id => $_ } map { $_->get_taxon } @{ $mrca->get_terminals };
				@tanglees = sort { $a cmp $b } grep { $_ ne $name } map { $_->get_name  } values %t;
			}
		}
		push @result, [ $name, $status, join(',',@tanglees), join(',',@ids), join(',',@meta) ];
	});
	return @result;
}

# makes a taxa block that contains all the distinct species, with their ID annotations
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
				$taxon_name = join ' ', @parts[0,1];
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

# reads extra metadata from spreadsheet, attaches records to taxa
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

# for each taxon, collects its stopnodes and indexes them with left and right counters
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
				
				# if two distinct taxa, it's a 'stopnode':
				# http://biophylo.blogspot.nl/2013/04/algorithm-for-distinguishing-polyphyly.html
				if ( scalar(keys(%merged)) >= 2 ) {
									
					# re-set the tally
					$node->set_generic( 'taxa' => {}, 'tanglees' => [ values %merged ] );				

					# tell all the subtended taxa about their stopnode				
					for my $taxon ( values %merged ) {
						my $stopnodes = $taxon->get_generic('stopnodes') || [];
						push @$stopnodes, $node;
						$taxon->set_generic( 'stopnodes' => $stopnodes );
					}
				}
				
				# continue the tally
				else {
					$node->set_generic( 'taxa' => \%merged );
				}
			}
		},
	);
}

# process command line / CGI arguments
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
		$verbosity  = INFO;
		$format     = $cgi->param('format');
		$separator  = $cgi->param('separator');
		$comments   = $cgi->param('comments');
		$whitespace = $cgi->param('whitespace');
		$trinomials = $cgi->param('trinomials');
		$infile     = $cgi->param('infile');
		$metadata   = $cgi->param('metadata');
		$quotes     = $cgi->param('quotes');
		$astsv      = $cgi->param('astsv');
		
		# TO DO: maybe make it so that users can also provide a location
		# for the input data, so that the monophylizer can be embedded
		# in BioVeL/Taverna workflows.
		$infh       = $cgi->upload('infile')->handle;
		$metafh     = $cgi->upload('metadata')->handle if $metadata;
		if ( CGI->can('Accept') ) {
			$astsv = ( CGI::Accept('text/plain') == 1.0 );
		}
		elsif ( CGI->can('CGI::accept') ) {
			$astsv = ( CGI::accept('text/plain') == 1.0 );			
		}
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

# writes an HTML table
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
