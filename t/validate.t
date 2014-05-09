#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use FindBin '$Bin';
use File::Spec;
use Data::Dumper;

# construct and verify required locations
my $script   = File::Spec->catfile( $Bin, '..', 'script', 'monophylizer.pl' );
my $datadir  = File::Spec->catdir( $Bin, '..', 'data' );
my $expected = File::Spec->catfile( $Bin, 'expected.tsv' );

ok( -x $script,   "is executable: $script" );
ok( -d $datadir,  "is a directory: $datadir" );
ok( -e $expected, "exists: $expected" );

# make command
my $command = $^X . ' ' . join( ' ', map { "-I$_" } @INC ) . ' ' . $script;

# read the file with expected output
my @exp;
{
	open my $fh, '<', $expected or die $!;
	@exp = parse_tsv($fh);
}

# read from the data directory
opendir my $dh, $datadir or die $!;
while( my $entry = readdir $dh ) {

	# verify files with nwk extension
	if ( $entry =~ /\.nwk$/ ) {
	
		# construct path for file
		my $infile = File::Spec->catfile( $Bin, '..', 'data', $entry );
		
		# run the script
		diag("going to run '$command -i $infile'");
		my $result = `$command -i $infile`;
		open my $fh, '<', \$result;
		
		# validate the output
		my @res = parse_tsv($fh);
		my @focal = grep { $_->{'Tree'} eq $entry } @exp;
		for my $f ( @focal ) {
			my $species = $f->{'Species'};
			my ($res) = grep { $_->{'Species'} eq $species } @res;
			my ( $obs, $exp ) = ( $res->{'Assessment'}, $f->{'Assessment'} );
			ok( $obs eq $exp, "Assessment for $species ($entry) - obs: $obs exp: $exp" );
			
			# validate whether the tanglees match
			if ( $f->{'Tanglees'} ) {
				my %expected = map { $_ => 1 } split /,/, $f->{'Tanglees'};
				my %observed = map { $_ => 1 } split /,/, $res->{'Tanglees'};
				is_deeply( \%observed, \%expected, "Tanglees" );
			}
		}		
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