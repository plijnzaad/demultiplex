#!/usr/bin/env perl
#
# Given a barcode file, demultiplexes a FASTQ file (on stdin) while potentially allowing for mismatches.
# For format of the barcode file, see testdata/testbarcodes.txt.
#
# To test, do e.g. 
#
#   ./demultiplex.pl -m 1 < testdata/one-mismatch.fastq  -b testdata/testbarcodes.txt -p DEMUL
# 
# written by <plijnzaad@gmail.com>

use strict;
use Getopt::Std;
use FileHandle;

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

my $barcodes_mixedcase = mismatch::readbarcodes_mixedcase($opt_b); ## eg. $h->{'AGCGtT') => 'M3'
my $barcodes = mismatch::mixedcase2upper($barcodes_mixedcase);     ## e.g. $h->{'AGCGTT') => 'M3'

my @mismatch_REs; $#mismatch_REs= ($allowed_mismatches);

for(my $i=1; $i<=$allowed_mismatches; $i++) { 
  my $re= mismatch::convert2mismatchREs(barcodes=>$barcodes_mixedcase, 
                                        allowed_mismatches =>$allowed_mismatches);# eg. $h->{'AGCGTT') =>  REGEXP(0x25a7788)
  $mismatch_REs[$i]=$re;
  ## note: leave $mismatch_REs[0] empty to avoid confusion (could use it for exact matches, but to slow)
}

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
RECORD:
while(1) { 
  my $record=<>;
  ### e.g.:  ^@NS500413:172:HVFHWBGXX:1:11101:4639:1062 1:N:0:CCGTCCAT$
  my ($foundcode)=(split(':', $record))[-1];
  $foundcode =~ s/[\n\r]*$//;
  $record .= <>; # sequence line
  $record .= <>; # '+'
  $record .= <>; # quality line
  
  my $lib;
 CASE:
  while(1) {
    $lib=$barcodes->{$foundcode};       # majority of cases
    if ($lib) {
      $nexact++;
      last CASE;
    }
    if (! $allowed_mismatches) {
      $nunknown++;
      $lib='UNKNOWN';
      last CASE;
    }
    my $correction;
    my $i;
  TRY:
    for($i=1; $i<= $allowed_mismatches; $i++) { 
      $correction=mismatch::rescue($foundcode, $mismatch_REs[$i]);
      last TRY if $correction;
    }
    if($correction) {
      $lib=$barcodes->{$correction};
      $nmismatched++;
      last CASE;
    } else { 
      $nunknown++;
      $lib='UNKNOWN';
      last CASE;
    }
    die "should not reach this point";
  }                                     # CASE
  $filehandles->{$lib}->print($record);
  last RECORD if (eof(STDIN) || !$record);
}                                       # RECORD
mismatch::close_outfiles($filehandles);

warn sprintf("exact: %s\nrescued: %s\nunknown: %s\n", 
             map { mismatch::commafy $_ } ($nexact, $nmismatched, $nunknown ));


