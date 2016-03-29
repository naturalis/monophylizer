#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO qw'parse unparse';

# process command line arguments
my ( $infile, $format, $out );
GetOptions(
	'infile=s' => \$infile,
	'format=s' => \$format,
	'out=s'    => \$out,
);

my $project = parse(
	'-format' => $format,
	'-file'   => $infile,	
	'-as_project' => 1,
);

print unparse(
	'-phylo'  => $project,
	'-format' => $out,
);