#!/usr/bin/perl
use strict;
use warnings;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

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

my $help;
GetOptions(
	'infile=s'    => \$infile,
	'verbose+'    => \$verbosity,
	'format=s'    => \$format,
	'separator=s' => \$separator,
	'comments'    => \$comments,
	'whitespace'  => \$whitespace,
	'trinomials'  => \$trinomials,
	'metadata=s'  => \$metadata,
	'help|?'      => \$help,
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
	'-keep_whitespace' => !!$whitespace,
);
$log->info("done reading tree: $tree");

# read the metadata spreadsheet, if provided
my %spreadsheet;
my $fieldcount;
if ( $metafh ) {
	$log->info("going to read additional metadata from $metadata");
	while(<$metafh>) {
		chomp;
		my @record = split /\t/, $_;
		$fieldcount = scalar @record;
		my $species = shift @record;
		$spreadsheet{$species} = \@record;
	}
	$log->info("done reading metadata");
}

# now traverse the tree
my %tipmeta;
my %seen;
my $counter = 0;
my $metacount;
$log->info("going to traverse the tree");
$tree->visit_depth_first(
	'-pre' => sub {
		my $node = shift;
		
		# get the name for each tip
		if ( $node->is_terminal ) {
			my $name = $node->get_name;
			
			# split the name into taxon name (first part) and subsequent metadata
			my ( $taxon, @meta ) = split /\Q$separator\E/, $name;
			if ( $taxon ) {

				# split trinomials
				my @parts = split /\s/, $taxon;
				if ( not $trinomials and scalar @parts > 2 ) {
					$log->info("will interpret '$taxon' as '$parts[0] $parts[1]'");
					$taxon = $parts[0] . ' ' . $parts[1];
				}
				elsif ( $trinomials and scalar @parts > 2 ) {
					$log->info("will retain subspecific epithet in '$taxon'");					
				}

				# store the metadata that was embedded in the tip label as
				# a record in a hash table
				$tipmeta{$taxon} = [] if not $tipmeta{$taxon};
				push @{ $tipmeta{$taxon} }, \@meta;
				$metacount = scalar @meta;
				$log->debug("found taxon '$taxon' with metadata '@meta'");

				# start building a mapping from species names to lists of tips,
				# once we get deeper in the tree we can then test whether these
				# lists are monophyletic
				$node->set_generic( 'tips' => { $taxon => [ $node ] } );
			}
			
			# somehow badly formatted name, emit warning
			else {
				$log->warn("couldn't process name '$name' with separator '$separator'");
			}
		}
		else {
			$node->set_generic( 'pre' => $counter++ );
		}
	},
	'-post' => sub {
		my $node = shift;
		
		# now process the deeper nodes
		if ( $node->is_internal ) {
			$node->set_generic( 'post' => $counter++ );
			my %tips;
			
			# merge the species => tips mappings of all children
			for my $child ( @{ $node->get_children } ) {
				my $taxa = $child->get_generic('tips');
				for my $taxon ( keys %{ $taxa } ) {
					if ( not $tips{$taxon} ) {
						$tips{$taxon} = $taxa->{$taxon};
					}
					else {
						push @{ $tips{$taxon} }, @{ $taxa->{$taxon} };
					}
				}
			}
			
			# more than 1 species subtended. in this situation the following
			# scenarios may apply:
			# 1. this is a deep node that subtends multiple, correctly lineage-sorted
			#    species
			# 2. this is a shallow node where an individual from species A pops up
			#    within a cluster of species B (i.e. species A is paraphyletic with
			#    respect to species B)
			# 3. this is a shallow node where species A and B are mangled together,
			#    (i.e. they are polyphyletic with respect to each other, so report
			#    all permutations)
			if ( 1 < scalar keys %tips ) {
				my $polynode;
				for my $taxon ( keys %tips ) {
					if ( $seen{$taxon} ) {
						$log->info("poly- or paraphyly detected for $taxon");
						push @{ $seen{$taxon} }, $node;
						$polynode++;
					}
					else {
						$seen{$taxon} = [ $node ];
					}
				}

				# reset the hash, we will evaluate later whether it was poly or para								
				$node->set_generic( 'tips' => {} );
				
				# store the tip labels if this is a polynode so that we can later
				# report the tangles
				$node->set_generic( 'tangles' => [ keys %tips ] );
			}
			
			# carry forward to deeper nodes
			else {
				$node->set_generic( 'tips' => \%tips );
			}
		}
	}	
);

# now do the final reporting
if ( $as_cgi_process ) {
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
		<h1>Results for $infile (click column headers to sort)</h1>
		<table class="sortable"><tr>
HEADER
	for my $i ( 0 .. ( $metacount + $fieldcount + 2 ) ) {
		print "<th>Column $i</th>";
	}
	print '</tr>';
}
for my $taxon ( sort { $a cmp $b } keys %tipmeta ) {
	my @row = ( $taxon );

	# now check to see if the taxon is mono/poly/paraphyletic
	my @ancestors = @{ $seen{$taxon} };
	if ( scalar @ancestors == 1 ) {
		$log->debug("'$taxon' appears to be monophyletic");
		push @row, 'monophyletic', '';
	}
	else {
		# I *think* we can now proceed as follows:
		# we first sort all the ancestors by their pre-order labeling. then, when we iterate
		# over them, the post-order label of the focal node should be larger than that of the
		# next one if the two are nested (i.e. paraphyletic), otherwise they are polyphyletic
		my @sorted = sort { $a->get_generic('pre') <=> $b->get_generic('pre') } @ancestors;
		my $jumps = 0;
		for my $i ( 0 .. $#sorted ) {
			if ( $i < $#sorted ) {
				my $pre      = $sorted[ $i     ]->get_generic('pre');
				my $nextpre  = $sorted[ $i + 1 ]->get_generic('pre');
				my $post     = $sorted[ $i     ]->get_generic('post');
				my $nextpost = $sorted[ $i + 1 ]->get_generic('post');
				if ( $pre + 1 != $nextpre && $post < $nextpost ) {
					$jumps++;
				}
			}
		}

		# i.e. the ancestors were not all nested
		if ( $jumps ) {
			$log->info("'$taxon' appears to be polyphyletic, seen $jumps jumps");
			push @row, 'polyphyletic';
		}
		else {
			$log->info("'$taxon' appears to be paraphyletic, seen $jumps jumps");
			push @row, 'paraphyletic';
		}

		# the sorted ancestors at this point also have the ancestor of the MRCA of the poly tips, we must ignore this
		my %tangles = map { $_ => 1 } grep { /\S/ } map { @{ $_->get_generic('tangles') } } @sorted[ 1 .. $#sorted ];
		delete $tangles{$taxon}; # delete self
		push @row, join ',', keys %tangles;
	}

	# first collapse all the species-level metadata
	my @records = @{ $tipmeta{ $taxon } };
	$log->debug("going to collapse ".scalar(@records)." metadata records");
	my $fieldcount = $#{ $records[0] };
	for my $i ( 0 .. $fieldcount ) {
		my %values;
		for my $record ( @records ) {
			$values{ $record->[$i] } = 1;
		}
		push @row, join ',', sort { $a cmp $b } keys %values;
	}

	# finally, append any pre-existing metadata from the spreadsheet
	push @row, @{ $spreadsheet{$taxon} } if $spreadsheet{$taxon};

	# print result
	if ( $as_cgi_process ) {
		s/\s/&nbsp;/g for @row;
		print "<tr><td>", join( "</td><td>", @row ), "</td></tr>\n";
	}
	else {
		print join( "\t", @row ), "\n";
	}
}

print "</table></body></html>" if $as_cgi_process;
