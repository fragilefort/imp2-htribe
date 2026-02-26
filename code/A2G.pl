#!/usr/bin/env perl

use strict;
use warnings;

my $fasta_file = $ARGV[0];
my $out = $ARGV[1];
my $base_index = $ARGV[2];

die "Usage: $0 <input.fa> <output.fa> <base_index[INT]>\n" unless @ARGV == 3;

open(my $fh, "<", $fasta_file) or die "Cannot open $fasta_file: $!";
open(my $oh, ">", $out) or die "Cannot open $out: $!";

my $seq = "";
my $header = "";


while(my $line = <$fh>) {
    $line =~ s/\R//g;
    if($line =~ /^>/) {
        $header = $line;
        next; 
    }
    $seq .= $line;
}

my $len_seq = length($seq);
print "Header is $header\n";
print "The length of the sequence is $len_seq\n" ;

my $base_toedit = substr($seq, $base_index, 1);
print "The base which will be edited is $base_toedit\n";

if ($base_toedit =~ /[Aa]/) {
    substr($seq, $base_index, 1) = "G";
} else {
    warn "The base you are trying to edit is not A\n";
}

$len_seq = length($seq);
my $edited_base = substr($seq, $base_index, 1);
print "The length of the sequence after editing is $len_seq\n";
print "The base after editing is $edited_base\n"; 
print $oh "$header\n";
print $oh "$seq\n";

close $fh;
close $oh;
