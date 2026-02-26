#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %option;
my $USAGE = "Usage: $0 -n normal_gene.fa -e edited_gene.fa -o outdir\n";
getopts('n:e:o:h', \%option);

die $USAGE if ($option{h} || !$option{n} || !$option{e} || !$option{o});

my $norm_fa = $option{n};
my $edit_fa = $option{e};
my $out_dir = $option{o};
my $read_len = 100;
my $fold = 1000; 
mkdir $out_dir unless -d $out_dir;

system("art_illumina -p -i $norm_fa -f $fold -l $read_len -m 200 -s 10 -o $out_dir/MASTER_NORM_");
system("art_illumina -p -i $edit_fa -f $fold -l $read_len -m 200 -s 10 -o $out_dir/MASTER_EDIT_");

my $total_reads = 10000;
my $norm_reads_count = $total_reads * 0.9;
my $edit_reads_count = $total_reads * 0.1;

my $norm_lines = $norm_reads_count * 4;
my $edit_lines = $edit_reads_count * 4;

for my $i (1..9) {
    my $out_r1 = "$out_dir/ctrl_DS${i}_R1.fq";
    my $out_r2 = "$out_dir/ctrl_DS${i}_R2.fq";

    if ($i <= 6) {
        system("head -n " . ($total_reads * 4) . " $out_dir/MASTER_NORM_1.fq > $out_r1");
        system("head -n " . ($total_reads * 4) . " $out_dir/MASTER_NORM_2.fq > $out_r2");
    } else {
        # 90% Normal + 10% Edited
        system("head -n $norm_lines $out_dir/MASTER_NORM_1.fq > $out_r1");
        system("head -n $norm_lines $out_dir/MASTER_NORM_2.fq > $out_r2");
        
        system("head -n $edit_lines $out_dir/MASTER_EDIT_1.fq >> $out_r1");
        system("head -n $edit_lines $out_dir/MASTER_EDIT_2.fq >> $out_r2");
    }
}
