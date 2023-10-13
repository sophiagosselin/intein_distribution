#!/usr/bin/perl -w
use strict;
use warnings;

#Opens the database of interest (in FASTA format), looks for the asc provided, then extracts those sequences.
#furthermore, this was designed to work on a proprietary database based on phagesDB phams
#So it will also take the sequences found, take the associated phams and ensure all pham members are also extracted
#returns a fasta file with all sequences of interest

#this was built originally to find: major_tail_protein
#ex ascession: >Makai_20_uniqueID_474 AU5 10050  major_tail_protein    474


#path to databse to search
my $database = $ARGV[0];
#ascession to look for
my $search_term = $ARGV[1];

#takes database file and finds all entries with the search term in ascession location
open(IN, "< $database");
open(OUT, "+> $search_term.fasta");
my $toggle = 0;
while(<IN>){
  chomp;
  if($_=~/\>/){
    if($_=~/$search_term/){
     $toggle=1;
     print OUT "$_\n";
    }
    else{
      $toggle=0;
      next;
    }
  }
  else{
    if($toggle == 1){
     print OUT "$_\n";
    }
    else{
     next;
    }
  }
}
close IN;
close OUT;
