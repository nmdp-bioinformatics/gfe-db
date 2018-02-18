#!/usr/bin/perl

use strict;
use warnings;
use Bio::AlignIO;

my $s_input_msf = shift @ARGV or die "NO INPUT FILE $! $0";
my $s_output_sth = shift @ARGV or die "NO OUTPUT FILE $! $0";

my $in = Bio::AlignIO->new(-file => $s_input_msf ,
                        -format => 'msf');

my $out = Bio::AlignIO->new(-file => ">".$s_output_sth,
                         -format => 'stockholm');

while ( my $aln = $in->next_aln ) {
    $out->write_aln($aln);
}

