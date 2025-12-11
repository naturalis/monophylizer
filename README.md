monophylizer
============
A web and command line tool to assess monophyly

Installation instructions
=========================
This tool depends on the Bio::Phylo toolkit, which needs to be installed from the 
comprehensive perl archive network here: http://search.cpan.org/dist/Bio-Phylo

Installation from CPAN is a standard system administration task for which ample 
documentation is available on the internet, and which is easy for UNIX-capable
system administrators and ICT support staffers. Consult them if you don't know 
how to do this.

To run the tool on the command line, do this:

	perl monophylizer.pl -infile newickFile.nwk > outputTable.tsv

Additional command line arguments (which you almost certainly won't need) are 
documented in the source code of monophylizer.pl

To run this tool as a web service, the following steps are necessary:

* Place monophylizer.pl in a location where it will be executed as a perl CGI script
by your web server. This will probably involve copying it into a folder called cgi-bin 
(linux/unix) or CGI-Executables (OSX). In addition, it may need to be set to "executable" 
and the first line of the script may need to be changed to point to the correct location
of the perl interpreter if it is not located at /usr/bin/perl

* Place monophylizer.html and sorttable.js in a location for plain text (not executable)
files.

* If the javascript file sorttable.js is not in the "document root", you will need to 
edit monophylizer.pl such that the output HTML can find it. This is at time of writing
at line 250.

Build status
============

Currently, the build status is:

[![CI](https://github.com/naturalis/monophylizer/actions/workflows/ci.yml/badge.svg)](https://github.com/naturalis/monophylizer/actions/workflows/ci.yml)

How to cite
===========
Pending publication of the method, this tool can be cited as:

**Marko, M., et al.**, 2016. Species-level para-and polyphyly in DNA barcode gene trees: strong operational bias in European Lepidoptera. 
_Systematic Biology_. **65**(6):1024-1040 doi:[10.1093/sysbio/syw044](https://doi.org/10.1093/sysbio/syw044)

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.16874.svg)](http://dx.doi.org/10.5281/zenodo.16874)

