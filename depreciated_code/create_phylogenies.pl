#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;

my $csv = $ARGV[0];
MAIN();

sub MAIN{
  my @clusters_of_interest = READIN_CSV($csv);
  mkdir("phylogenies");

  foreach my $coi (@clusters_of_interest){
    #move phylogenies into seperate directories
    mkdir("phylogenies\/$coi");
    copy("blast_searches\/$coi\_extein.fasta","phylogenies\/$coi\/$coi\_extein.fasta");
    copy("clusters\/$coi","phylogenies\/$coi\/$coi\_intein.fasta");

    #align sequences
    my $extein_alignment = MUSCLE("phylogenies\/$coi\/$coi\_extein.fasta");
    my $intein_alignment = MUSCLE("phylogenies\/$coi\/$coi\_intein.fasta");

    #create phylogenies
    IQTREE($extein_alignment);
    IQTREE($intein_alignment);
  }
}

sub IQTREE{
  #takes aligned fasta file
  #returns phylogenetic tree

  my $alignment = shift;
  system("iqtree2 -s $alignment -m MFP -B 1000 -T AUTO");

}

sub MUSCLE{
  #takes protein fasta as input
  #returns muscle aligned fasta file

  my $infasta = shift;
  system("muscle -align $infasta -output $infasta.aligned");

  return("$infasta.aligned");
}

sub READIN_CSV{
  #readin csv of clusters of interest
  #return array of those clusters

  my $csv_file=shift;
  my @csv_content;

  open(my $IN, "< $csv_file");
  while(<$IN>){
    chomp;
    push(@csv_content,$_);
  }
  close $IN;

  return(@csv_content);
}
