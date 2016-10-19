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

## Using the code

To start using it, do 

  git clone https://github.com/plijnzaad/demultiplex.git

This will create a directory `demultiplex`. Make sure that directory is on your `$PATH`. 
If you use the `mismatch.pm` library, make sure that directory is part of your `$PERLLIB`.
To get a 'usage message', run the script with a `-h` flag.

This code used to reside in https://plijnzaad@github.com/plijnzaad/phtools, subdirectory `ngs/demultiplex`

I'm interested in feedback and improvements, please report issues and/or fork it and send me pull requests.
