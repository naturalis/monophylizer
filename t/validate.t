#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use FindBin '$Bin';
use File::Spec;

# construct and verify required locations
my $script   = File::Spec->catfile( $Bin, '..', 'script', 'monophylizer.pl' );
my $datadir  = File::Spec->catdir( $Bin, '..', 'data' );
my $expected = File::Spec->catfile( $Bin, 'expected.tsv' );

ok( -x $script, "is executable: $script" );
ok( -d $datadir, "is a directory: $datadir" );
ok( -e $expected, "exists: $expected" );

# read the file with expected output
my %exp;
{
	open my $fh, '<', $expected or die $!;
	%exp = parse_tsv($fh);
}

# read from the data directory
opendir my $dh, $datadir or die $!;
while( my $entry = readdir $dh ) {

	# verify files with nwk extension
	if ( $entry =~ /\.nwk$/ ) {
	
		# construct path for file
		my $infile = File::Spec->catfile( $Bin, '..', 'data', $entry );
		
		# run the script
		my $result = `$script -i $infile`;
		open my $fh, '<', \$result;
		my %res = parse_tsv($fh);
	}
}

sub parse_tsv {
	my $fh = shift;
	my @header;
	my @result;
	while(<$fh>) {
		chomp;
		my @record = split /\t/, $_;
		if ( not @header ) {
			@header = @record;
		}
		else {
			my %record;
			for my $i ( 0 .. $#header ) {
				$record{$header[$i]} = $record[$i];
			}
			push @result, \%record;
		}
	}
	return @result;
}