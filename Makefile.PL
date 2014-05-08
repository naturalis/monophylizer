#!/usr/bin/perl
use strict;
use warnings;
use ExtUtils::MakeMaker;

my %parms = (
    'NAME'         => 'Bio::Monophylizer', # EU::MM apparently now wants Package::Name
    'AUTHOR'       => 'Rutger Vos',
    'PL_FILES'     => {},
    'EXE_FILES'    => [ 'script/monophylizer.pl' ],
    'VERSION_FROM' => 'script/monophylizer.pl',
    'PREREQ_PM'    => { 'Bio::Phylo' => 0 },
    'LICENSE'      => 'perl',
    'ABSTRACT'     => 'Script to assess monophyly in BOLD trees',
    'clean'        => {},    
    'dist'         => {
        'COMPRESS' => 'gzip -9f',
        'SUFFIX'   => 'gz',
        'TARFLAGS' => '-c -v -f',
    },
);

if ( $ExtUtils::MakeMaker::VERSION ge '6.46' ) {
    $parms{'META_MERGE'} = {
        'resources' => {
            'homepage'    => 'http://biophylo.blogspot.com/',
            'bugtracker'  => 'https://github.com/naturalis/monophylizer/issues',
            'repository'  => 'git://github.com/naturalis/monophylizer.git',
            'license'     => 'http://dev.perl.org/licenses/',
            'MailingList' => 'mailto:bio-phylo@googlegroups.com',
        },
    };
}

WriteMakefile( %parms );