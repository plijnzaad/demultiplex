#!/usr/bin/env perl
# See usage

use strict;
use Getopt::Std;
use FileHandle;

use mismatch;

use vars qw($opt_h $opt_b $opt_m $opt_p $opt_o $opt_e $opt_g);

my $version=mismatch::getversion($0);
my @args=@ARGV;

$opt_e=  '^cbc=([ACTGN]+)';             # select the barcode 

my $Usage="Usage:

   samtools view bigfile.bam | $0 -b barcodes.txt [ -g barcodegroups.txt ] [ -e $opt_e ] [-m mismatches] [ -p outputprefix ] [ -o outputdir ] 

Given a barcode file, demultiplexes a SAM file (on stdin) where the
barcode is expected in the QNAME (first field), and extracted from it
using a parenthesized regular expression indicated by the -e option (the
original use case is to split single-cell data by their cell bar
codes). As with the demultiplex-fastq.pl script, mismatches can be allowed
using the -m option.

Output files are written in *BAM* format, with all the original headers. 

For format of the barcode file, see testdata/testbarcodes.txt.

If a file barcodegroup is given (option -g), all barcode id's are 
replaced with their group ids, allowing the merging of several barcodes
into  larger groups. The format of the group file is simply 
  barcode_id \\s  group_id
As with the barcode file, comments (/#.*$/) and empty lines are ignored.

To test, do e.g. 

  (see demultiplex.pl \@FIX)

NOTE: the script does *not* check if mismatched barcodes are unambiguous!
Use edit-distance.pl and/or edit-distance-matrix.pl for that. To fix
ambiguous barcodes of a badly picked set of barcodes, run
check-barcodes.pl (it can adjust the barcodes so that mismatches in
ambiguous positions are disallowed).

written by <plijnzaad\@gmail.com> ";

if ( !getopts("b:p:o:m:g:h") || ! $opt_b ||  $opt_h ) {
    die $Usage; 
}

die "-m option missing " unless defined($opt_m);

warn "Running $0, version $version, with args @args\n";

my  $allowed_mismatches = $opt_m;

my $barcodes_mixedcase = mismatch::readbarcodes_mixedcase($opt_b); ## eg. $h->{'AGCGtT') => 'M3'
my $barcodes = mismatch::mixedcase2upper($barcodes_mixedcase);     ## e.g. $h->{'AGCGTT') => 'M3'

my $mismatch_REs=[]; $#{$mismatch_REs}= ($allowed_mismatches);

for(my $i=1; $i<=$allowed_mismatches; $i++) { 
  my $re= mismatch::convert2mismatchREs(barcodes=>$barcodes_mixedcase, 
                                        allowed_mismatches =>$allowed_mismatches);# eg. $h->{'AGCGTT') =>  REGEXP(0x25a7788)
  $mismatch_REs->[$i]=$re;
  ## note: leave $mismatch_REs->[0] empty to avoid confusion (could use it for exact matches, but too slow)
}

$barcodes_mixedcase=undef;

my $groups=undef;
$groups=mismatch::read_groups($opt_g) if $opt_g;

my @files=();


if ($groups) { 
  my $occur; map {  $occur->{$_}++ } values %$groups;
  @files= (keys %$occur, 'UNKNOWN');
}else { 
  @files= (values %$barcodes, 'UNKNOWN');
}

my $filehandles=mismatch::open_outfiles(outdir=>$opt_o, 
                                        prefix=>$opt_p,
                                        type=>'bam',
                                        files=>[@files]);      # opens M3.bam, UNKNOWN.bam etc.

my $nexact=0;
my $nmismatched=0;                         # having at most $mismatch mismatches
my $nunknown=0;

my $barcode_re = qr/$opt_e/;

my $nrefseqs=0;
my $warned=0;

my $stdin=FileHandle->new_from_fd(0, "<") || die "stdin: $!";

## lastly, process the actual input:
my $results=mismatch::demultiplex(type=>'sam',
                               input=>$stdin, 
                               outputs=>$filehandles,
                               barcodes=>$barcodes, 
                               mismatch_REs=>$mismatch_REs,
                               groups=>$groups,
                               barcode_re=>$barcode_re);

mismatch::close_outfiles($filehandles);

warn sprintf("exact: %s\nrescued: %s\nunknown: %s\n", 
             map { mismatch::commafy $_ } (map {$results->{$_}} qw(nexact nmismatched nunknown)));
