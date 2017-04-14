#!/usr/bin/env perl
=head1 NAME

    expected-observed-gfe.pl

=head1 SYNOPSIS
    

=head1 AUTHOR     Mike Halagan <mhalagan@nmdp.org>
    
    Bioinformatics Scientist
    3001 Broadway Stree NE
    Minneapolis, MN 55413
    ext. 8225

=head1 DESCRIPTION


=head1 CAVEATS

  *** This is not production ready and has not been extensively tested ***

=head1 LICENSE

    Copyright (c) 2017 National Marrow Donor Program (NMDP)

    This library is free software; you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation; either version 3 of the License, or (at
    your option) any later version.

    This library is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.
 
    You should have received a copy of the GNU Lesser General Public License
    along with this library;  if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.

    > http://www.gnu.org/licenses/lgpl.html

=head1 VERSIONS
    

=head1 TODO
    

=head1 SUBROUTINES

=cut
use strict;    # always
use warnings;  # or else
use FindBin;
use Data::Dumper;
use FindBin;
use Math::Round;
use GFE::Client;
use GFE::Annotate;
use XML::Simple;

my $s_xml_file = shift @ARGV or die "No XML file provided!";
my $s_outname  = shift @ARGV or die "No output name provided!";
my $s_oldname  = shift @ARGV or die "No output name provided!";

my %h_done;
if(defined $s_oldname){
  print STDERR "Loading..\n";
  open(my $fh,"<",$s_oldname) or die "CANT OPEN FILE $! $0";
  while(<$fh>){
    chomp;
    my($s_allele,@x) = split(/,/,$_);
    $h_done{$s_allele}++;
  }
  close $fh;
}
my $n = scalar keys %h_done;
print STDERR "$n alleles already finished\n";
my $verbose = 1;

my $o_gfe = GFE->new();
my $o_structures = $o_gfe->structures;
my $o_annotate   = $o_gfe->annotate;
my $o_client = $o_gfe->client;

my $s_expected = $s_outname.".expected.txt";
my $s_observed = $s_outname.".observed.txt";


# 
# download hla.xml 
# $url = "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip";
#
my $file = `cat $s_xml_file`;

print STDERR "parsing..." if $verbose;
my $alleles = XMLin($file);
print STDERR "done\n" if $verbose;

my %GF;  # Gene Feature
foreach my $allele (keys %{$alleles->{"allele"}}) {
  next if defined $h_done{$allele};
  # $allele='HLA-A*01:01:08'
  my $a = $alleles->{allele}->{$allele};
  my $s = $a->{sequence};
  my $nucsequence= $s->{nucsequence};
  foreach my $feature (keys %{$s->{feature}}) {
    # $feature = "Exon 8";
    next unless $feature =~/(Exon|Intron|UTR)/;
    my $start = $s->{feature}->{$feature}->{SequenceCoordinates}->{start};
    my $end = $s->{feature}->{$feature}->{SequenceCoordinates}->{end};
    my $nseq = substr($nucsequence, $start-1, $end-$start+1);
    $GF{$allele}{$feature} = $nseq;
  }
  $GF{$allele}{sequence} = $nucsequence;
}

my $num_alleles = scalar keys %GF;

my $nums = 1;
my %h_skipped;
open(my $fh_expected,">",$s_expected) or die "CANT OPEN FILE $! $0";
foreach my $allele (sort keys %GF) {
  $nums++;
  my %h_accesion = ();
  my %h_seqs = ();
  my $locus = (split /\*/, $allele)[0];
  foreach my $feature (keys %{$GF{$allele}}) {
    my $term;
    my $rank;
    if ($feature=~ /(\D+) (\d+)/ )  {
      # Exon 2
      # Intron 3
      $term = lc $1; # lower case
      $rank = $2;
    } elsif ($feature=~/3' UTR/) {
      $term = "three_prime_UTR";
      $rank = 8;
    } elsif ($feature=~/5' UTR/) {
      $term = "five_prime_UTR";
      $rank = 1;
    } else {
      #warn "unknown feature: $feature";
      next;
    }
    # format subject to change; 
    # print join (' ', $locus, $term, $rank, $GF{$allele}{$feature}), "\n";
    $h_seqs{$term}{$rank} = $GF{$allele}{$feature};
    my $n_accesion = $o_client->getAccesion($locus,$term,$rank,$GF{$allele}{$feature});
    $h_accesion{$term}{$rank} = $n_accesion;
  }

  if(!defined $o_annotate->order->{$locus}){
    print STDERR "Not valid $locus : Skipped $allele\n";
    $h_skipped{$allele}++;
    next;
  }

  if(scalar keys %h_accesion == 0){
    print STDERR "Skipped $allele\n";
    $h_skipped{$allele}++;
    next;
  }

  my @a_gfe;
  my @a_terms;
  push(@a_terms,join(";","Sequence",$GF{$allele}{sequence}));
  foreach my $term_rank (@{$o_structures->getStruct($locus)}) {
      my ($term, $rank) = split /:/, $term_rank;
      $h_accesion{$term}{$rank} = defined $h_accesion{$term}{$rank} ? $h_accesion{$term}{$rank} : '0';
      my $s_seq = defined $h_seqs{$term}{$rank} ? $h_seqs{$term}{$rank} : '';
      push(@a_terms,join(";",$term_rank,$h_accesion{$term}{$rank},$s_seq ));
      push(@a_gfe, $h_accesion{$term}{$rank}); 
  }

  if($nums % 100 == 0){
    my $per = $nums / $num_alleles;
    print STDERR $per."% done ($nums / $num_alleles)\n";
  }

  ## Create GFE
  my $s_gfe        = join('w',$locus, join('-', @a_gfe));
  print $fh_expected join(",",$allele,$s_gfe,@a_terms),"\n";

}
close $s_expected;

my $nums = 1;
my %h_skipped;
open(my $fh_observed,">",$s_observed) or die "CANT OPEN FILE $! $0";
foreach my $allele (sort keys %GF) {
  my $locus = (split /\*/, $allele)[0];
  my $rh_gfe = $o_gfe->getGfe($locus,$GF{$allele}{sequence});

  if(!defined $o_annotate->order->{$locus}){
    print STDERR "Not valid $locus : Skipped $allele\n";
    $h_skipped{$allele}++;
    $nums++;
    next;
  }

  if(defined $$rh_gfe{Error}){
    print STDERR "Skipped $allele ".$$rh_gfe{Error}{Message}."\n";
    $h_skipped{$allele}++;
    $nums++;
    next;
  }

  my %h_acces;
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

  if($nums % 100 == 0){
    my $per = $nums / $num_alleles;
    print STDERR $per."% done ($nums / $num_alleles)\n";
  }

  print $fh_observed join(",",$allele,$$rh_gfe{gfe},@a_terms),"\n";
  $nums++;
}
close $fh_observed;


exit 0;