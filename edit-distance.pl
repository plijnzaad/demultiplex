#!/usr/bin/env perl
use strict;
use mismatch;
use Getopt::Std;
## use Text::Levenshtein qw(distance);     

use vars qw($opt_h $opt_m);

my $Usage="Find near-duplicates among barcodes. 

For list of strings on stdin (format: name TAB barcode), print all the
sets of strings with distance less than or equal NUMBER to each other.
Mismatched letters are 'highlighted' in lowercase.

Usage: 
  $0 -m NUMBER  [<] barcodefile.txt > output
";

if ( !getopts("hm:") || $opt_h ) {
    die $Usage; 
}

die "-m option missing " unless defined($opt_m);

my $limit= $opt_m || 2;


my  $strings={};

### read codes
LINE:
while(<>) { 
  s/#.*//;
  s/[ \t\n\r]*$//g;
  next LINE unless $_;
  my ($name, $code)=split(' ');            # e.g. 'G7 \t CCAACAAT'
  if (! $code) {                           # just one column of strings
    $code=$name;
    $name= "(line $.)";
  }
  die "duplicate string: '$code' (named '$name', line $.,)" if $strings->{$code};
  $strings->{$code}=$name;
}                                       # while

my @strings = sort keys %$strings;


print "# distance with edit distance <= $limit:\n";
SEQ:
for(my $i=0; $i<@strings; $i++) { 
  my @s=@strings;
  my $s=$s[$i];
  splice(@s, $i, 1);
  ## my @d = distance($s , @s); ### Levenshtein edit distance, but we don't allow indels
  my @d = map { mismatch::hammingdist($_, $s); } @s;
  my @hits = grep( $d[$_] <= $limit, 0..$#d);
  next SEQ unless @hits;
  print $strings->{$s}. "\t$s\n"; 
  for(my $h=0; $h<@hits; $h++) {
    my $hs=$s[ $hits[$h] ];
    my $mm=mismatch::format_mm($s, $hs);
    print "$strings->{$hs}\t$mm (d=$d[ $hits[$h] ])\n";
  }
  print "\n";
}
