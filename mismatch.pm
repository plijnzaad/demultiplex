
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
  ## expects a barcode ID and sequence (3-12 nt long) on a single line separated by a single tab
  my ($file)=@_;
  my $barcodeids={};
  my $barcodes_mixedcase = {};
  my $uppercase_codes={};
  my $len=undef;
  my $nlowercase=0;

  open(FILE, "$file") or die "Barcode '$file': $!";
LINE:
  while(<FILE>) {
    chomp;
    s/#.*//;
    next LINE unless $_;
    my ($barcodeid, $code)=split("\t",$_);  # e.g. 'G7 \t CCAACAAT'
    my $l=length($code);

    if ( defined($len) && $l != $len) { 
      die "Barcode $code has length $l, expected uniform length $len, file $file line $.";
    } else {
      $len=$l;
    }
    my $regexp='^[ACGT]{3,12}$';
    die "Barcode $code  does not match /$regexp/i, file $file, line $." unless $code =~ /$regexp/i;
    # (3,12 are a wild guess at sanity)

    $nlowercase += ($code =~ /[a-z]/) ;

    die "Barcode id '$barcodeid' not unique" if $barcodeids->{$barcodeid}++;
    die "Barcode '$code' not unique" if $uppercase_codes->{"\U$code"}++;
    $barcodes_mixedcase->{$code}=$barcodeid;
  }                                     # while LINE
  close(FILE);
  warn "Found $nlowercase barcodes containing lower case letters. These will be uppercased and will not be allowed to mismatch at the lowercased positions" 
      if $nlowercase;
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
## takes hash with barcodes (e.g. $h->{'AGCGtT') => 'M3' ) and the allowed number of mismatches
## returns per barcode the mismatch regular expresson e.g. $h->{'AGCGTT') => REGEXP(0x25a7788)
## The resulting map only contains uppercase barcodes, as this is needed for mapping it to the
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
  my ($barcodes, $max_mismatches)= map {$args->{$_}} qw(barcodes max_mismatches);

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
  ##mark the fixed positons with !
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
  for(my $i=1; $i<=$max_mm; $i++) { 
    ## set up all possible combinations of mismatch positions (usually just 1, since max_mm usually 1)
    ## combine (Math::Combinatorics) returns all unique unordered combinations of $i mismatch positions in a set of length($code) barcode positions 
    my @pos_sets = combine($i, 0..$#code);
  COMB:
    foreach my $pos_set ( @pos_sets ) { 
      ## replace the mismatch positions with '.' (regexp for "any character") using splicing (yay)
      my @mm=@code;
      @mm[ @$pos_set ] = split(//, '.' x int(@$pos_set) ); 
      my $mm_re=join("", @mm);
      for my $i (0 .. $#fixed) { 
        if ($fixed[$i] eq '!' && $mm[$i] eq '.') { 
          ## warn "regexp $mm_re conflicts with a fixed position, skipped\n";
          next COMB;
        }
      }
      push(@mmcodes, $mm_re);
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
  
  my($nexact, $nunknown, $nrescued, $statsperbarcode, $statspermm); 
  my($nrefseqs, $warned);               # only used for bam

  foreach my $code (keys %$barcodes )  { 
    $statsperbarcode->{$code}=[];
    $statspermm->{$code}={};
  }
  my $filehandles=$outputs;

RECORD:
  while(1) { 
    my $foundcode;                      
    my $record=<$input>;
    ## ($foundcode and $record are the only two variables needed)
    if ($type eq 'fastq') { 
      ### e.g.:  ^@NS500413:172:HVFHWBGXX:1:11101:4639:1062 1:N:0:CCGTCCAT$
      ### this code will NOT work for FASTQ file produced by Casava < 1.8 see https://en.wikipedia.org/wiki/FASTQ_format
      ### descriptions used in FASTQ files from NCBI/EBI will end up in $foundcode 
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
        $statsperbarcode->{$foundcode}->[0]++;
        last CASE;
      }
      if (! $mismatch_REs) {
        $nunknown++;
        $lib='UNKNOWN';
        last CASE;
      }
      my $correction;
      my $nmismatches;
    TRY:
      for($nmismatches=1; $nmismatches < @$mismatch_REs; $nmismatches++) { 
        $correction=mismatch::rescue($foundcode, $mismatch_REs->[ $nmismatches ]);
        last TRY if $correction;
      }
      if($correction) {
        $lib=$barcodes->{$correction};
        $nrescued->[$nmismatches]++;
        $statsperbarcode->{$correction}->[$nmismatches]++;
        $statspermm->{$correction}->{$foundcode}++;
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
  {nexact=>$nexact, nrescued=>$nrescued, nunknown=>$nunknown, 
   statsperbarcode=>$statsperbarcode, statspermm=>$statspermm};
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
  #expects per line a space separated barcode Id and group
  my($file)=@_;
  open(FILE, $file) || die "$0: $file: $!";
  my $groups={};

  while(<FILE>) { 
    s/#.*//;
    #s/[\r\n]*$//;
    chomp;
    next unless /\S+\s+\S+/;
    my($barcode,$group)=split(' ',$_);
    die "barcode $barcode not unique in group file $file, line $.," if $groups->{$barcode};
    $groups->{$barcode}=$group;
  }
  close(FILE);
  $groups;
}                                       # sub read_group

sub print_statsperbarcode { 
  my $args = ref $_[0] eq 'HASH' ? shift : {@_}; # args: file, stats, max_mismatches, barcodes
  my ($stats, $mismatches, $max_mismatches, $barcodes)= 
      map {$args->{$_}} qw(stats mismatches max_mismatches barcodes);

  ## overall stats:
  my $file="counts-overall.txt";
  warn "Creating file $file ...\n";
  open(OUT, "> $file")  || die "$file: $!";
  my $plain="exact\t" . join("\t", map { "${_}mm"; } 1..$max_mismatches);
  my $perc="%exact\t" . join("\t", map { "%".$_."mm"; } 1..$max_mismatches);
  print OUT "#id\tcode\t$plain\t$perc\n";
 CODE:
  foreach my $code (sort keys %$barcodes) { 
    my $id=$barcodes->{$code};
    my $nexact=$stats->{$code}->[0] || 0;
    print OUT "$id\t$code\t$nexact\t"; 
    my $total=$nexact;

  I:
    for(my $i=1; $i<=$max_mismatches; $i++) { 
      my $n=$stats->{$code}->[$i] || 0;
      print OUT mismatch::commafy($n) . "\t";
      $total += $n;
    }
    ## now percentages:
    printf OUT "%.1f\t", 100*$nexact/$total;
  I:
    for(my $i=1; $i<=$max_mismatches; $i++) { 
      my $n=$stats->{$code}->[$i] || 0;
      printf OUT "%4.1f\t", 100*$n/$total . "\t";
    }
    print OUT "\n";
  }                                     # CODE
  close(OUT);

  ### now per mismatch
  $file="counts-permismatch.txt";
  warn "Creating file $file ...\n";
  open(OUT, "> $file")  || die "$file: $!";
  print OUT "#id\tbarcode\tmismatched barcodes with counts\n";
CODE:
  foreach my $code (sort keys %$barcodes) { 
    my $id=$barcodes->{$code};
    print OUT "$id\t$code\t";
    my $mms=$mismatches->{$code};
    for my $m (sort keys %$mms) {
      print OUT format_mm($code, $m) . ": $mms->{$m}\t";
    }
    print OUT "\n";
  }                                     # CODE
  close(OUT);
}                                       # sub print_statsperbarcode
1;
