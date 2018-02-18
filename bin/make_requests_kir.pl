#!/usr/bin/env perl
##############################################################################
# SCRIPT NAME:	mk_requests_kir
# DESCRIPTION:	extract features from KIR.dat
#
# DATE WRITTEN: 2016-02-12
# WRITTEN BY:   Martin Maiers
#
##############################################################################
use strict;    # always
use warnings;  # or else
use FindBin;
use Math::Round;
use lib "$FindBin::Bin/../lib";
use GFE;
use GFE::Client;
use GFE::Annotate;
use XML::Simple;

my $verbose = 1;

my $o_gfe = GFE->new();
my $o_structures = $o_gfe->structures;
my $o_annotate = $o_gfe->annotate;
my $o_client = $o_gfe->client;


#
# donload ENA formatted KIR.dat from IPD kir
#
#my $file = "./KIR.dat.txt";
#open FILE, $file or die "$!: $file";


my %K;
my $ac = "";
my $de = "";
my $state = "";
my $start = "";
my $stop = "";
my $intron;
my $exon;
my $UTR;
my $seq;
while(<>) {
  chomp;
  
  next unless /^( |AC|DE|FT)/;
  if (/^ /) {
    $K{$ac}{seq}.= $1 if /([actg ]+)/;
    $K{$ac}{seq}=~s/ //g;

  }
   
  if (/^AC   (\S+);/) {
    $seq = "";
    $ac = $1; 
    $intron=1; 
    $exon=1; 
    $UTR=1; 
    $start=1; 
    $stop=1; 
  }
  if (/^DE   (\S+),/) {
    $K{$ac}{de} = $1;
  }

  if ($state eq "intron") {
    $intron = $1 if /(\d+)/;
    $state = "";
    $K{$ac}{"intron$intron"}{start} = $start;
    $K{$ac}{"intron$intron"}{stop}   = $stop;
  }
  if ($state eq "exon") {
    $exon = $1 if /(\d+)/;
    $state = "";
    $K{$ac}{"exon$exon"}{start} = $start;
    $K{$ac}{"exon$exon"}{stop}   = $stop;
  }

  if(/^FT   UTR\s+(\d+)\.\.(\d+)/) {
    $start = $1;
    $stop  = $2;
    $K{$ac}{"UTR$UTR"}{start}  = $start;
    $K{$ac}{"UTR$UTR"}{stop}   = $stop;
    $UTR++;
  }
  if(/^FT   exon\s+(\d+)\.\.(\d+)/) {
    $state = "exon";
    $start = $1;
    $stop  = $2;
  }
  if(/^FT   intron\s+(\d+)\.\.(\d+)/) {
    $state = "intron";
    $start = $1;
    $stop  = $2;
  }
}

use Data::Dumper;
my @a = keys %K;
# print join("\n",@a),"\n";
# #print Dumper(%K),"\n";
# print STDERR "Finished loading!\n";

# my %h_skipped;
# my $s_expected = "observed.kir.txt";
# open(my $fh_expected,">",$s_expected) or die "CANT OPEN FILE $! $0";
# my $nums = 0;
# foreach my $ac (sort keys %K) {

#   $nums++;
#   my %h_accesion = ();
#   my %h_seqs = ();
#   my $locus = (split /\*/, $K{$ac}{de})[0];
#   foreach my $term (sort keys %{$K{$ac}}) {
#     next unless $term =~/^(exon|intron|UTR)/;

#     my $seq = $K{$ac}{seq};
#     my $start = $K{$ac}{$term}{start};
#     my $stop = $K{$ac}{$term}{stop};

#     #print $seq, "\n";
#     my $t;
#     my $r;
#     if ($term=~/(\D+)(\d+)/) {
#       $t= $1;
#       $r= $2;
#       if ($t eq "UTR" && $r eq "1") {

#         # need to address the case where there is only a 3' UTR
#         if ($start eq "1") {
#           $t = 'five_prime_UTR';
#         } else {
#           $t = 'three_prime_UTR';
#         }
#       }
#       if ($t eq "UTR" && $r eq "2") {
#         $t = 'three_prime_UTR';
#       }
#     }
#     my $sequence = uc substr($seq, $start, ($stop-$start+1));

#     $h_seqs{$t}{$r} = $sequence;
#     my $n_accesion = $o_client->getAccesion($locus,$t,$r,$sequence);
#     $h_accesion{$t}{$r} = $n_accesion;

#     if(!defined $o_annotate->order->{$locus}){
#       print STDERR "Not valid $locus : Skipped $K{$ac}{de}\n";
#       $h_skipped{$K{$ac}{de}}++;
#       next;
#     }

#     if(scalar keys %h_accesion == 0){
#       print STDERR "Skipped $K{$ac}{de}\n";
#       $h_skipped{$K{$ac}{de}}++;
#       next;
#     }
#   }

#   if(!defined $o_annotate->order->{$locus}){
#     print STDERR "Not valid $locus : Skipped $K{$ac}{de}\n";
#     $h_skipped{$K{$ac}{de}}++;
#     next;
#   }
#   next if scalar (keys %h_accesion) == 0;

#   my @a_gfe;
#   my @a_terms;
#   push(@a_terms,join(";","Sequence",$K{$ac}{seq}));
#   foreach my $term_rank (@{$o_structures->getStruct($locus)}) {
#       my ($term, $rank) = split /:/, $term_rank;
#       $h_accesion{$term}{$rank} = defined $h_accesion{$term}{$rank} ? $h_accesion{$term}{$rank} : '0';
#       my $s_seq = defined $h_seqs{$term}{$rank} ? $h_seqs{$term}{$rank} : '';
#       push(@a_terms,join(";",$term_rank,$h_accesion{$term}{$rank},$s_seq ));
#       push(@a_gfe, $h_accesion{$term}{$rank}); 
#   }

#   my $s_gfe        = join('w',$locus, join('-', @a_gfe));
#   #print STDERR join(",",$K{$ac}{de},$s_gfe,@a_terms),"\n";
#   print $fh_expected join(",",$K{$ac}{de},$s_gfe,@a_terms),"\n";
# }
# close $s_expected;

my %h_skipped;
my $s_expected = "observed.kir.txt";
open(my $fh_expected,">",$s_expected) or die "CANT OPEN FILE $! $0";
my $nums = 0;
foreach my $ac (sort keys %K) {

  $nums++;
  my %h_acces = ();
  my $allele = $K{$ac}{de};
  my $locus = (split /\*/, $K{$ac}{de})[0];

  if(!defined $o_annotate->order->{$locus}){
    print STDERR "Not valid $locus : Skipped $K{$ac}{de}\n";
    $h_skipped{$K{$ac}{de}}++;
    next;
  }

  my $rh_gfe = $o_gfe->getGfe($locus,$K{$ac}{seq});
  if(defined $$rh_gfe{Error}){
    print STDERR "Skipped $allele ".$$rh_gfe{Error}{Message}."\n";
    $nums++;
    next;
  }

  my %h_seq;
  foreach my $rh_structure (@{$$rh_gfe{structure}}){
      $h_acces{$$rh_structure{term}}{$$rh_structure{rank}} = $$rh_structure{accession};
      $h_seq{$$rh_structure{term}}{$$rh_structure{rank}}   = $$rh_structure{sequence};
  }

  my @a_terms;
  my $f_aligned = $$rh_gfe{aligned};
  push(@a_terms,join(";","Sequence",$$rh_gfe{fullgene}{sequence}));
  foreach my $term_rank (@{$o_structures->getStruct($locus)}) {
      my ($term, $rank) = split /:/, $term_rank;
      $h_acces{$term}{$rank} = defined $h_acces{$term}{$rank} ? $h_acces{$term}{$rank} : '0';
      my $s_seq = defined $h_seq{$term}{$rank} ? $h_seq{$term}{$rank} : '';
      push(@a_terms,join(";",$term_rank,$h_acces{$term}{$rank},$s_seq ));
  }

  print $fh_expected join(",",$allele,$$rh_gfe{gfe},@a_terms),"\n";
  $nums++;

}
close $s_expected;



exit 0;

