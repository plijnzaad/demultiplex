
package mismatch;

### Usage: see demultiplex-{fastq,sam}.pl

use strict;

use Math::Combinatorics;
use Regexp::Optimizer;
use File::Basename;
use FileHandle;

sub readbarcodes {
  ## utility function to read barcodes, returns hash with eg. $barcodes->{'AGCGtT') => 'M3' }
  ## Note that lower case letters (used for disallowing specific mismatches) are
  ## still there (and won't match actual barcodes).
  my ($file)=@_;
  my $barcodeids={};
  my $barcodes_mixedcase = {};
  my $uppercase_codes={};
  my $len=undef;

  open(FILE, "$file") or die "Barcode '$file': $!";
LINE:
  while(<FILE>) {
    s/[\n\r]*$//g;
    s/#.*//;
    next LINE unless $_;
    my ($barcodeid, $code)=split(' ');  # e.g. 'G7 \t CCAACAAT'
    my $l=length($code);

    if ( defined($len) && $l != $len) { 
      die "Barcode $code has length $l, expected uniform length $len, file $file line $.";
    } else {
      $len=$l;
    }

    die "Barcode $code  does not match /^[ACGT]{3,12}$/i, file $file, line $." unless $code =~ /^[ACGT]{3,12}$/i;
    # (3,12 are a wild guess at sanity)

    warn "barcode $barcodeid contains lower case letters, these will be uppercased and will not be allowed to mismatch" 
        if $code =~ /[a-z]/ ;

    die "Barcode id '$barcodeid' not unique" if $barcodeids->{$barcodeid}++;
    die "Barcode '$code' not unique" if $uppercase_codes->{"\U$code"}++;
    $barcodes_mixedcase->{$code}=$barcodeid;
  }                                     # while LINE
  close(FILE);
  $barcodes_mixedcase;
}                                       # readbarcodes

sub mixedcase2upper { 
  ## utility function to convert the mixed case hash (which is used for the mismatch regular expressions) to an uppercased
  ## hash
  my ($mixed) = @_;
  my $barcodes={};
  for my $code ( keys %$mixed ) { $barcodes->{"\U$code"}=$mixed->{$code}}
  $barcodes;
}

sub convert2mismatchREs {
## takes hash with barcodes (e.g. $h->{'AGCGtT') => 'M3' ) and returns
## e.g. $h->{'AGCGTT') => REGEXP(0x25a7788). The resulting map only
## contains uppercase barcodes, as this is needed for mapping it to the
## output file.  The hash returned contains, per barcode, one regexp
## representing all possible mismatches of that barcode.  In the values
## (i.e. regexps), lowercase letters (if any) are uppercased and the
## regexp does not allow these letters to mismatch.

  my $args = ref $_[0] eq 'HASH' ? shift : {@_}; # args: barcodes, allowed_mismatches
  my $o=Regexp::Optimizer->new;

  my $mm_REs={};
  for my $code (keys %{$args->{barcodes}}) {
    my @res= _getmismatch_REs($code, $args->{allowed_mismatches}); # empty if allowed_mismatches==0
    my $r='^'.join("|", @res).'$';
    $r=$o->optimize(qr/$r/);
    $mm_REs->{"\U$code"}= $r;           # just one big regexp. Note the uppercasing
  }                                     # for $code
  $mm_REs;  
}                                       # convert2mismatchREs

sub get_mismatchREs {
  ### set up array of regexp for increasing numbers of mismatches, to be
  ### tested in turn.  They are to be tested in turn. Element 0
  ### (corresponding to 'no mismatches') is deliberately left undefined
  ### to avoid confusion

  my $args = ref $_[0] eq 'HASH' ? shift : {@_};
  my ($barcodes, $max_mismatches)= map {$args->{$_}} qw(barcode max_mismatches);

  return undef if ($max_mismatches==0);

  my $mismatch_REs=[]; 
  $#{$mismatch_REs}= ($max_mismatches);

  for(my $i=1; $i<=$max_mismatches; $i++) { 
    my $re= mismatch::convert2mismatchREs(barcodes=>$barcodes, 
                                          allowed_mismatches =>$i);
    ## eg. $h->{'AGCGTT') =>  REGEXP(0x25a7788)
    $mismatch_REs->[$i]=$re;
  }
  $mismatch_REs;
}                                       # get_mismatchREs


sub rescue { 
  ### return the barcode without mismatches (not its ID!)
  my($foundcode, $mm_REs)=@_;

  foreach my $code (keys %$mm_REs) {
    my $re=$mm_REs->{$code};
    return $code if $foundcode =~ $re;
  }
  return undef;
}                                       # rescue

sub safe_rescue { 
  my($foundcode, $mm_REs)=@_;

  my @found=();
  foreach my $code (keys %$mm_REs) {
    my $re=$mm_REs->{$code};
    push(@found,   $code) if $foundcode =~ $re;
  }
  return @found;
}                                       # rescue


sub format_mm { 
  ## given an real barcode and the found one (with mismatches), make all the
  ## mismatched letters lowercase
  my($orig, $mm)=@_;
  my @orig=split('', $orig);
  my @mm=split('', $mm);
  my @new;
  for(my $i=0; $i<@orig; $i++) {
    push(@new, ($orig[$i] eq  $mm[$i])? $mm[$i] : "\L$mm[$i]");
  }
  join("", @new);
}


sub _getmismatch_REs {
  ## for one barcode, set up the regular expressions that allows mismatches
  my($code, $max_mm)=@_;

  return () if ! $max_mm;

  my @fixed = ();
  if ($code =~ /[a-z]/)  {
    my $fixed= $code;
    $fixed =~ s/[a-z]/!/g;
    @fixed = split(//, $fixed);
    $code = "\U$code";
  }

  my @mmcodes=();
  my(@code)=split(//, $code);

  ## set up array of arrays with '.' where to do the replacements:
  for(my $i=0; $i<$max_mm; $i++) { 
    my @combs = combine(($i+1), 0..$#code);
  COMB:
    foreach my $comb ( @combs ) { 
      my @mm=@code;
      @mm[ @$comb ] = split(//, '.' x int(@$comb) ); # yay, splicing
      my $mm=join("", @mm);
      for my $i (0 .. $#fixed) { 
        if ($fixed[$i] eq '!' && $mm[$i] eq '.') { 
          warn "regexp $mm conflicts with a fixed position, skipped\n";
          next COMB;
        }
      }
      push(@mmcodes, $mm);
    }
  }
  @mmcodes;
}                                       # _getmismatch_REs

sub hammingdist {                       
## honestly stolen from http://www.perlmonks.org/?node_id=500244
  length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] );
}

sub getversion {
  my($path)=@_;
  my ($fullpath)=`which $path`;
  my ($script,$dir) = fileparse($fullpath);
  my $version=`cd $dir && git describe --match 'v[0-9]*' --tags --dirty`;
  chomp($version);
  $version='UNKNOWN' unless $version;
  $version;
}                                       # getversion

sub commafy {
  # insert comma's to separate powers of 1000
  my($i)=@_;
  my $r = join('',reverse(split('',$i)));
  $r =~ s/(\d{3})/$1,/g;
  $r =~ s/,$//;
  join('',reverse(split('',$r)));
}

sub demultiplex {
  my $args = ref $_[0] eq 'HASH' ? shift : {@_};
  my ($type, $input, $outputs, $barcodes, $mismatch_REs, $groups, $barcode_re)=
      map {$args->{$_}} qw(type input outputs barcodes mismatch_REs groups barcode_re);

  die "unknown type '$type', must be fastq or bam" if ($type ne 'fastq' && $type ne 'bam');
  
  my($nexact, $nunknown, $nrescued); 
  my($nrefseqs, $warned);               # only used for bam

  my $filehandles=$outputs;

RECORD:
  while(1) { 
    my $foundcode;                      
    my $record=<$input>;
    ## ($foundcode and $record are the only two variables needed)
    if ($type eq 'fastq') { 
      ### e.g.:  ^@NS500413:172:HVFHWBGXX:1:11101:4639:1062 1:N:0:CCGTCCAT$
      $foundcode=(split(':', $record))[-1];
      $foundcode =~ s/[\n\r]*$//;
      $record .= <$input>; # sequence line
      $record .= <$input>; # '+'
      $record .= <$input>; # quality line
    } else { 
      ### sam file, header line:
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
      ## else: sam file, read line:
## e.g. ^NS500413:188:H3M3WBGXY:1:11101:10124:1906:cbc=TACCTGTC:umi=TTCGAC \t 0 \t GLUL__chr1 \t 3255 \t 25 \t 76M \t 
      my($qname,$flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen,
         $seq, $qual, @optionals)=split("\t", $record);
      
      for my $part (split(":", $qname)) {
        $foundcode=$1 if $part =~ $barcode_re;
      }
      die "could not find barcode in QNAME '$qname', expected /$barcode_re/, line $." unless $foundcode;
    }
    ### at this point we need and have just $foundcode and $record
    my $lib;
  CASE:
    while(1) {
      $lib=$barcodes->{$foundcode};       # majority of cases
      if ($lib) {
        $nexact++;
        last CASE;
      }
      if (! $mismatch_REs) {
        $nunknown++;
        $lib='UNKNOWN';
        last CASE;
      }
      my $correction;
      my $i;
    TRY:
      for($i=1; $i < @$mismatch_REs; $i++) { 
        $correction=mismatch::rescue($foundcode, $mismatch_REs->[$i]);
        last TRY if $correction;
      }
      if($correction) {
        $lib=$barcodes->{$correction};
        $nrescued++;
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
    last RECORD if ( $input->eof() || !$record );
  }                                       # RECORD
  {nexact=>$nexact, nrescued=>$nrescued, nunknown=>$nunknown};
}                                         # sub demultiplex

sub open_infile {
  die "not used nor tested";
  my($file)=@_;
  my $fh=FileHandle->new();
  if ($file =~ /\.gz/) { 
    $fh->open("zcat $file | ", "r")  or die "'$file': $!";
  } else { 
    $fh->open("< $file")  or die "'$file': $!";
  }
  $fh;
}

sub close_infile {
  die "not used nor tested";
}

sub open_outfiles { 
  my $args = ref $_[0] eq 'HASH' ? shift : {@_};
  my ($outdir, $prefix, $type, $files)=map {$args->{$_}} qw(outdir prefix type files);
  my(@libs)=@$files;
  my $fhs={};

  die "Output directory $outdir: $!" if ($outdir && !(-d $outdir && -w $outdir));

  for my $lib (@libs) { 
    my $name;
    my $fh;
    if ($type =~ /fastq/) {
      $name=sprintf("%s.fastq.gz", $lib);
      $name="$prefix$name" if $prefix;
      $name="$outdir/$name" if $outdir;
      $fh = FileHandle->new("| gzip -n > $name");
    } elsif ($type eq 'bam' ) { 
      $name=sprintf("%s.bam", $lib);
      $name="$prefix$name" if $prefix;
      $name="$outdir/$name" if $outdir;
      $fh = FileHandle->new(" | samtools view - -h -b > $name");
    } else { 
      die "open_outfiles: unknown type '$type' requested";
    }
    die "library $lib, file $name: $!" unless $fh;
    warn "Creating/overwriting file $name ...\n";
    $fhs->{$lib}=$fh;
  }
  $fhs;
}                                       # open_outfiles

sub close_outfiles {
  my($fhs)=@_;
  for my $lib (keys %$fhs) {
    $fhs->{$lib}->close() or die "could not close (or open?) demultiplexed bam file for library $lib; investigate";
  }
}

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
}                                       # sub read_group


1;
