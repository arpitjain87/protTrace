#!/usr/bin/perl -w -I /home/CIBIV/ingo/perlmod
use strict;
use Getopt::Long;

# PROGRAMNAME: degapper.pl

# AUTHOR: INGO EBERSBERGER, ingo.ebersberger@univie.ac.at

# PROGRAM DESCRIPTION:

# DATE: Thu May 29 11:11:16 CEST 2008


# DATE LAST MODIFIED: 08.07.2013
# 	mdofied the output name
######################## start main #############################
my $help;
my $infile;
my $outfile;
my %final;
my $limit = 0.5;
my @sum = qw();
GetOptions ("h" => \$help,
	    "in=s" => \$infile,
	    "limit=s" => \$limit,
	    "out=s" => \$outfile);
if ($help or !defined $infile) {
	die "degapper.pl -in=<> [-limit=<>] [-out=<>] [-h]

DESCRIPTION:
degapper.pl processes an alignment by removing all alignment columns where less than 'limit' species are represented by a valid amino acid.
OPTIONS:
-in: provide the filename of the alignment in fasta format
-out: provide the outfile name (optional)
-limit: determines the minimum fraction of species that have to be represented by a amino acid in an alignment column. DEFAULT: 0.5";
}

if (!defined $outfile) {
	$outfile=$infile;
	$outfile=~s/(.*)\.(.*)/$1.proc.$2/;
	print "$1\n";
	print "outfile is $outfile\n";
}

my @seqs = `less $infile |sed -e 's/>//'`;
chomp @seqs;
my %smkeep = @seqs;
my %sm = @seqs;
my $speccount = scalar(keys %sm);
## translate all aminoacids into a 1, all other characters into a 0
## and sum up
for (keys %sm) {
    $sm{$_} =~ s/[^FLIMVSPTAYHQNKDECWRG]/0/ig;
    $sm{$_} =~ s/[A-Z]/1/gi;
    my @loc = split //, $sm{$_};
    for (my $i = 0; $i < @loc; $i++) {
	$sum[$i] += $loc[$i];
    }
}

## extract the alignment
for (my $i = 0; $i < @sum; $i++) {
    if ($sum[$i]/$speccount > $limit) {
	for (keys %sm) {
	    $final{$_} .= substr($smkeep{$_}, $i, 1);
	}
    }
}
open (OUT, ">$outfile") or die "could not open outfile $outfile\n";
for (my $i = 0; $i < @seqs; $i+= 2) {
    print OUT ">$seqs[$i]\n$final{$seqs[$i]}\n";
}
close OUT;

