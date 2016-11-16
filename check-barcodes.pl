#!/usr/bin/env perl
# See usage

use strict;
use Getopt::Std;
use FileHandle;
use File::Basename;

use mismatch;
use Algorithm::Combinatorics;

use vars qw($opt_h $opt_m);

my $Usage="Usage:

   $0 barcodes.txt -m NMISMATCHES

written by <plijnzaad\@gmail.com>
";

if ( !getopts("m:h") || ! $opt_m || $opt_h || ! @ARGV) {
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

warn "Will now check all possible 5 ^ $len i.e. " . 5 ** $len . " barcodes\n";

my $iter= Algorithm::Combinatorics::tuples_with_repetition([ qw(A C G T N)] , $len);

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
  
  if (@codes>1) {
    print "$w: ". join(' ', @codes) . "\n";
  }

  for my $b (@codes) { 
    push @{$badcodes->{$w}}, @codes;
  }
}

warn "stats:\n";
warn "unknown: $nunknown\n";
warn "unique: $nunique\n";
warn "ambiguous: $nambiguous\n";
warn "done\n";
1;
