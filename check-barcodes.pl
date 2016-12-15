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
that have ambiguous potential barcode misreads.  The output (to stdout) is 

MISREAD : BARCODE1 BARCODE2 [ BARCODE3 etc.]  

where BARCODE1 etc. are the real barcodes, with letters mismatching the
MISREAD 'highlighted' in lowercase. Statistics are written to stderr.

When the -o FILE option is used, it will calculate new barcodes that are
written to FILE. These new barcodes are unambiguous while allowing for
mismatches. These barcodes have lowercase letters in places where
mistmatches cannot be tolerated.  This output can be used by the
demultiplex-{sam,fastq}.pl scripts (they use the same convention:
mismatches are never allowed on lowercase letters).  Note that running the
current script on this new file will therefore not find any ambiguity
anymore. For a testcase, see the file testdata/mm1-ambiguous-codes.txt

written by <plijnzaad\@gmail.com>
";

if ( !getopts("m:o:h") || $opt_h  || @ARGV!= 1) {
    die $Usage; 
}
my  $allowed_mismatches = $opt_m;

my $barcodefile=$ARGV[0];

my $barcodes_mixedcase = mismatch::readbarcodes_mixedcase($barcodefile); ## eg. $h->{'AGCGtT') => 'M3'

my $barcodes = mismatch::mixedcase2upper($barcodes_mixedcase);     ## e.g. $h->{'AGCGTT') => 'M3'
my $mismatch_REs = mismatch::convert2mismatchREs(barcodes=>$barcodes_mixedcase, 
                                                 allowed_mismatches =>$allowed_mismatches);# eg. $h->{'AGCGTT') =>  REGEXP(0x25a7788)

# $barcodes_mixedcase=undef;
my $id2code={};
## invert the table for lookup later on
for my $code (keys %$barcodes_mixedcase) { 
  $id2code->{  $barcodes_mixedcase->{$code} }=$code;
}

my $len=length( (keys %$barcodes)[0]  );

warn "Will now check all possible 4 ^ $len i.e. " . 4 ** $len . " barcodes for ambiguous misreads,
allowing for $allowed_mismatches mismatches\n";

my $iter= Algorithm::Combinatorics::tuples_with_repetition([ qw(A C G T) ] , $len);

my $nunknown=0;
my $nexact=0;
my $nunique=0;
my $nambiguous=0;
my $badcodes={};

sub merge_codes { 
  ## do something like qw(AcGTT ACgTT ACgTt ) => "AcgTt"
  my (@mm)=@_;

  return $mm[0]  if(@mm==1);

  my @w= map { [ split('', $_) ] } @mm;

  my $len=int(@{$w[0]});

  my @final=();

  for(my $i=0; $i<$len; $i++) { 
    my $h={};
    my @col=map { $_->[$i];} @w ;
    for my $col  ( @col ) { $h->{$col}++; }
    my $n= int(keys %$h);
    push(@final,  ($n==1) ? $col[0] : "\L$col[0]");
  }
  join("", @final);
}

print "# misreads that map ambiguously to mismatched barcodes:\n";
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

##  for my $code (sort keys %$badcodes) { 
  for my $id (sort keys %$id2code ) { 
    my $code=$id2code->{$id};           # mixed case!
    if(exists ($badcodes->{$code}))  {
      my @mm=keys %{$badcodes->{$code}};
      print FILE "$id\t". merge_codes(@mm) . "\t # was: $code\n";
    } else {
      print FILE "$id\t$code\n";
    }
  }
  close(FILE);
}


warn "\nstats:\n";
warn "no mismatch: $nexact\n";
warn "mismatched but unique: $nunique\n";
warn "mismatched ambiguous: $nambiguous\n";
warn "unknown: $nunknown\n";
warn "done\n";
1;
