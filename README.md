# demultiplex

Code for demultiplexing fastq files, and a library for flexibly allowing mismatches during barcode recognition.

## Contents

 * demultiplex.pl - the actual demultiplexer. Input is read from stdin, output is written to `*.fastq.gz`
 * edit-distance.pl - for spotting pairs of barcodes that are too close  in terms of Hamming distance
 * edit-distance-matrix.pl - print all the pairwise Hamming distance (to help choose or optimize barcodes)
 * mismatch.pm - library that deals with the mismatching
 * testdata - test data

## Mismatching

If the barcodes that were used are distinct enough, you can recover some
of the barcodes that are lost due to sequencing error (the default is to
allow 1 mismatch).  This is safe as long as the minimum edit distance
amongst all the barcodes is 2 or higher. If most of the distances are >=
2 but some are only 1, you can still allow most mismatches in these
barcodes, apart from the positions inside the barcodes where these pairs
of barcodes do not differ enough. This is done by specifying these
positions using *lowercase* letters. In contrast to the uppercase
letters, these are never allowed to mismatch. (Any lowercase letters in
the barcode file is uppercased before the matching; the lowercase is
only used to conveniently specifying the fixed position).

Note that if _k_ mistmatches are to be allowed, the minimum edit
distance among all barcodes must be at least _k+2_ (since otherwise the
mismatched code is exactly halfway between two exact barcodes).

The script barcode-check.pl can be used to check for clashed between
barcodes with mismatches.

## Using the code

To start using it, do 

`  git clone https://github.com/plijnzaad/demultiplex.git`

This will create a directory `demultiplex`. Make sure that directory is
on your `$PATH`.  If you use the `mismatch.pm` library, make sure that
directory is part of your `$PERLLIB`. The easiest way is to copy the
`env.sh.example` script to e.g. `env.sh`, adjust its contents, and
source (using `source env.sh` or `. env.sh`) it prior to running the
scripts or using the library. For a testrun, do e.g.

`
$ mkdir testoutput
$ ./demultiplex-fastq.pl -m 2 -b testdata/testbarcodes.txt -o testoutput  < testdata/nomismatches.fastq 

`
To get a usage message on a script, run the script with an `-h` flag.

This code used to reside in https://plijnzaad@github.com/plijnzaad/phtools, subdirectory `ngs/demultiplex`

I'm interested in feedback and improvements, please report issues and/or fork it and send me pull requests.
