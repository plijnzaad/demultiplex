#!/usr/bin/env perl
# See usage

use strict;
use Getopt::Std;
use FileHandle;
use File::Basename;

use mismatch;
use Algorithm::Combinatorics;

use vars qw($opt_h $opt_m $opt_o);

$opt_m=1;

my $Usage="Usage:

   $0  [ -m NMISMATCHES ] [-o uniqified-barcodes.txt ]   barcodes.txt

Given a barcodefile (format: id \\t barcode \\n), find all the barcodes
that have ambiguous potential barcode misreads.  The output is 

MISREAD : BARCODE1 BARCODE2 [ BARCODE3 etc.]  

where BARCODE1 etc. are the real barcodes, with letters mismatching the
MISREAD 'highlighted' in lowercase.

written by <plijnzaad\@gmail.com>
";

if ( !getopts("m:o:h") || $opt_h) {
    die $Usage; 
}
my  $allowed_mismatches = $opt_m;

my $barcodefile=$ARGV[0];

my $barcodes_mixedcase = mismatch::readbarcodes_mixedcase($barcodefile); ## eg. $h->{'AGCGtT') => 'M3'

my $barcodes = mismatch::mixedcase2upper($barcodes_mixedcase);     ## e.g. $h->{'AGCGTT') => 'M3'
my $mismatch_REs = mismatch::convert2mismatchREs(barcodes=>$barcodes_mixedcase, 
                                                 allowed_mismatches =>$allowed_mismatches);# eg. $h->{'AGCGTT') =>  REGEXP(0x25a7788)

$barcodes_mixedcase=undef;

my $len=length( (keys %$barcodes)[0]  );

warn "Will now check all possible 4 ^ $len i.e. " . 4 ** $len . " barcodes for ambiguous misreads,
allowing for $allowed_mismatches mismatches\n";

my $iter= Algorithm::Combinatorics::tuples_with_repetition([ qw(A C G T) ] , $len);

my $nunknown=0;
my $nexact=0;
my $nunique=0;
my $nambiguous=0;
my $badcodes={};

WORD:
while(my $inst=$iter->next()) {
  my $w=join("",@$inst);                # potential mismatched
  
  if (exists( $barcodes->{$w})) { 
    $nexact++;
    next WORD;
  }

  my @codes=mismatch::safe_rescue($w, $mismatch_REs);

  $nunknown +=  (@codes==0);
  $nunique +=  (@codes==1);
  $nambiguous +=  (@codes>1);
  
  if (@codes>1) {               # ambiguous
    my @mm = map { mismatch::format_mm($w, $_); } @codes;
    print "$w: ". join(' ', @mm) . "\n";
    
    for my $code (@codes) { 
      my $mm=mismatch::format_mm($w, $code);
      $badcodes->{$code}{$mm}++;
    }
  }
}                                       # WORD

if ($opt_o) { 
  open(FILE, "> $opt_o") || die "$opt_o: $!";

  die "not yet ready";
  close(FILE);
}


warn "\nstats:\n";
warn "no mismatch: $nexact\n";
warn "mismatched but unique: $nunique\n";
warn "mismatched ambiguous: $nambiguous\n";
warn "unknown: $nunknown\n";
warn "done\n";
1;
