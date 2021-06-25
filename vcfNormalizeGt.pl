#!/usr/bin/env perl

use warnings;
use strict;

my @sortedToInIx;

while (<>) {
  chomp; chomp;
  if (/^#CHROM/) {
    my @cols = split("\t");
    my %gtInIx;
    for (my $i = 9;  $i < @cols;  $i++) {
      $gtInIx{$cols[$i]} = $i;
    }
    my @sortedSamples = sort @cols[9..$#cols];
    for (my $i = 0;  $i < @sortedSamples;  $i++) {
      $sortedToInIx[$i] = $gtInIx{$sortedSamples[$i]};
    }
    print STDERR "Got " . scalar(@sortedSamples) . " samples, first $sortedSamples[0], "
      . "last $sortedSamples[-1]\n";
    print join("\t", "#POS", "REF", "ALT", @sortedSamples) . "\n";
  } elsif (! /^#/) {
    my @w = split("\t");
    # For comparison we only want pos, ref, and sorted list of alt alleles.
    my ($pos, $ref, $altStr) = ($w[1], $w[3], $w[4]);
    my @alts = split(",", $altStr);
    my $sortedAltStr = join(",", sort @alts);
    print join("\t", $pos, $ref, $sortedAltStr);
    # Output genotypes in sorted-sample, translated from numbers to base values
    for (my $i = 0;  $i < @sortedToInIx;  $i++) {
      my $ixIn = $sortedToInIx[$i];
      my $gt = $w[$ixIn];
      if (! defined $gt) {
        die "No gt for output sample $i, input column $ixIn";
      }
      my $allele = $gt;
      if ($gt eq "0") {
        $allele = $ref;
      } elsif ($gt =~ m/^\d+$/) {
        $allele = $alts[$gt-1];
        if (! defined $allele) {
          die "Bad gt '$gt' for alts $altStr at pos $pos, output sample $i, input column $ixIn";
        }
      }
      print "\t$allele";
    }
    print "\n";
  }
}
