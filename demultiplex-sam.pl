#!/usr/bin/env perl
# See usage

use strict;
use Getopt::Std;
use FileHandle;

use mismatch;

use vars qw($opt_h $opt_b $opt_m $opt_p $opt_o $opt_e $opt_g);

my $version=mismatch::getversion($0);
my @args=@ARGV;

$opt_e=  '^cbc=([ACTGN]+)';

my $Usage="Usage:

   samtools view bigfile.bam | $0 -b barcodes.txt [ -g barcodegroups.txt ] [ -e $opt_e ] [-m mismatches] [ -p outputprefix ] [ -o outputdir ] 

Given a barcode file, demultiplexes a SAM file (on stdin) where the
barcode is expected in the QNAME (first field), and extracted from it
using a parenthesized regular expression indicated by the -e option.
As with the original demultiplex script, mismatches can be allowed
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

Original copied from demultiplex.pl, 3bc1490195 (2016-10-19 17:40:41)
written by <plijnzaad\@gmail.com> ";

if ( !getopts("b:p:o:m:g:h") || ! $opt_b ||  $opt_h ) {
    die $Usage; 
}

die "-m option missing " unless defined($opt_m);

warn "Running $0, version $version, with args @args\n";

my  $allowed_mismatches = $opt_m;

my $barcodes_mixedcase = mismatch::readbarcodes_mixedcase($opt_b); ## eg. $h->{'AGCGtT') => 'M3'
my $barcodes = mismatch::mixedcase2upper($barcodes_mixedcase);     ## e.g. $h->{'AGCGTT') => 'M3'
my $mismatch_REs = mismatch::convert2mismatchREs(barcodes=>$barcodes_mixedcase, 
                                                 allowed_mismatches =>$allowed_mismatches);# eg. $h->{'AGCGTT') =>  REGEXP(0x25a7788)
$barcodes_mixedcase=undef;

my $groups=undef;
$groups=read_groups($opt_g) if $opt_g;

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

## lastly, process the actual input:
RECORD:
while(1) { 
  my $record=<>;

  if ($record =~ /^@/) {                # header line, needed by all files
    for my $lib (keys %$filehandles) { 
      $filehandles->{$lib}->print($record);
    }
    $nrefseqs += ($record =~ /^\@SQ/);
    next RECORD;
  }
  ### @@@FIX: at this point we should insert add a @PG record to the bam headers ...

  if ( $nrefseqs ==0 && !$warned++ ) {
    warn "*** expected to find reference sequences in the sam headers (the \@SQ records)\n";
    warn "*** be sure to use output from samtools -h\n";
  }

## e.g. ^NS500413:188:H3M3WBGXY:1:11101:10124:1906:cbc=TACCTGTC:umi=TTCGAC \t 0 \t GLUL__chr1 \t 3255 \t 25 \t 76M \t 
  my($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen,
     $seq, $qual, @optionals)=split("\t", $record);

  my $foundcode;
  for my $part (split(":", $qname)) {
    $foundcode=$1 if $part =~ $barcode_re;
  }
  die "could not find barcode in QNAME '$qname', expected /$barcode_re/, line $." unless $foundcode;
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
    my $correction = mismatch::rescue($foundcode, $mismatch_REs);
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
  $lib= $groups->{$lib} if $groups;
  $lib = 'UNKNOWN' unless $lib;

  $filehandles->{$lib}->print($record);
  last RECORD if (eof(STDIN) || !$record);
}                                       # RECORD
mismatch::close_outfiles($filehandles);

warn sprintf("exact: %s\nrescued: %s\nunknown: %s\n", 
             map { mismatch::commafy $_ } ($nexact, $nmismatched, $nunknown ));


sub read_groups { 
  #return hash mapping barcode to group
  my($file)=@_;
  open(FILE, $file) || die "$0: $file: $!";
  my $groups={};

  while(<FILE>) { 
    s/#.*//;
    s/[\r\n]*$//;
    next unless /\S+\s+\S+/;
    my($barcode,$group)=split(' ',$_);
    die "barcode $barcode not unique in group file $file, line $.," if $groups->{$barcode};
    $groups->{$barcode}=$group;
  }
  close(FILE);
  $groups;
}
