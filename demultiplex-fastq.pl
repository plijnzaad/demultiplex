#!/usr/bin/env perl
#
# Given a barcode file, demultiplexes a FASTQ file (on stdin) while potentially allowing for mismatches.
# For format of the barcode file, see testdata/testbarcodes.txt.
#
# To test, do e.g. 
#
#   ./demultiplex-fastq.pl -m 1 < testdata/one-mismatch.fastq  -b testdata/testbarcodes.txt -p DEMUL
# 
# written by <plijnzaad@gmail.com>

use strict;
use Getopt::Std;
use FileHandle;
## use IO::File;

use mismatch;

use vars qw($opt_h $opt_b $opt_m $opt_p $opt_o);

my $version=mismatch::getversion($0);
my @args=@ARGV;


my $Usage="Usage:

   zcat bigfile.fastq.gz | $0 -b barcodes.txt -m mismatches [ -p outputprefix ] [ -o outputdir ] 

NOTE: the script does *not* check if mismatched barcodes are unambiguous!
Use edit-distance.pl and/or edit-distance-matrix.pl for that. To fix
ambiguous barcodes of a badly picked set of barcodes, run
check-barcodes.pl (it can adjust the barcodes so that mismatches in
ambiguous positions are disallowed).


";

if ( !getopts("b:p:o:m:h") || ! $opt_b ||  $opt_h ) {
    die $Usage; 
}

die "-m option missing " unless defined($opt_m);

warn "Running $0, version $version, with args @args\n";

my  $allowed_mismatches = $opt_m;

my $barcodes_mixedcase = mismatch::readbarcodes($opt_b); ## eg. $h->{'AGCGtT') => 'M3'
my $barcodes = mismatch::mixedcase2upper($barcodes_mixedcase);     ## e.g. $h->{'AGCGTT') => 'M3'

my $mismatch_REs=mismatch::get_mismatchREs(barcodes=>$barcodes_mixedcase, max_mismatches=>$allowed_mismatches);
$barcodes_mixedcase=undef;

my @files=(values %$barcodes, 'UNKNOWN');

my $filehandles=mismatch::open_outfiles(outdir=>$opt_o, 
                                        prefix=>$opt_p,
                                        type=>'fastq',
                                        files=>[@files]);      # opens M3.fastq.gz, UNKNOWN.fastq.gz, etc.

my $nexact=0;
my $nmismatched=0;                         # having at most $mismatch mismatches
my $nunknown=0;

## lastly, process the actual input:

my $stdin=FileHandle->new_from_fd(0, "<") || die "stdin: $!";

my $results=mismatch::demultiplex(type=>'fastq',
                                  input=>$stdin, 
                                  outputs=>$filehandles,
                                  barcodes=>$barcodes, 
                                  mismatch_REs=>$mismatch_REs);

mismatch::close_outfiles($filehandles);

warn sprintf("exact: %s\n", mismatch::commafy($results->{nexact}));
for(my $i=1; $i<=$#{$results->{nrescued}}; $i++) {
  warn sprintf("%d mismatches: %s\n", $i, mismatch::commafy($results->{nrescued}->[$i]));
}
warn sprintf("unknown: %s\n", mismatch::commafy($results->{nunknown}));
