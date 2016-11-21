#!/usr/bin/env perl
# See usage

use strict;
use Getopt::Std;
use FileHandle;
use File::Basename;

use mismatch;
use Algorithm::Combinatorics;

use vars qw($opt_h $opt_m);

$opt_m=1;

my $Usage="Usage:

   $0  [ -m NMISMATCHES ]  barcodes.txt

Given a barcodefile (format: id \\t barcode \\n), find all the barcodes
that have ambiguous potential barcode misreads.  The output is 

MISREAD : BARCODE1 BARCODE2 [ BARCODE3 etc.]  

where BARCODE1 etc. are the real barcodes, with letters mismatching the
MISREAD 'highlighted' in lowercase.

written by <plijnzaad\@gmail.com>
";

if ( !getopts("m:h") || $opt_h) {
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
my $nunique=0;
my $nambiguous=0;
my $badcodes={};

while(my $inst=$iter->next()) {
  my $w=join("",@$inst);
  
  my @codes=mismatch::safe_rescue($w, $mismatch_REs);

  $nunknown +=  (@codes==0);
  $nunique +=  (@codes==1);
  $nambiguous +=  (@codes>1);
  
  if (@codes>1) {               # ambiguous
    my @mm = map { mismatch::format_mm($w, $_); } @codes;
    print "$w: ". join(' ', @mm) . "\n";
  }

  for my $b (@codes) { 
    push @{$badcodes->{$w}}, @codes;
  }
}

warn "\nstats:\n";
warn "unknown: $nunknown\n";
warn "unique: $nunique\n";
warn "ambiguous: $nambiguous\n";
warn "done\n";
1;
