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

my @phams;
#takes database file and finds all entries with the search term in ascession location
open(IN, "< $database");
my $toggle = 0;
while(<IN>){
  chomp;
  if($_=~/\>/){
    if($_=~/$search_term/){
      my($pham)=($_=~/\>.*\_UID\_\d+\s+.*?\s+(\d+?)\s+\w+\s+\d+$/);
      push(@phams,$pham);
      $toggle=1;
    }
    else{
      $toggle=0;
      next;
    }
  }
  else{
    next;
  }
}
close IN;

#takes list of phams and ensures these are extracted as well.
open(OUT, "+> $search_term.fasta");
open(IN, "< $database");
while(<IN>){
  chomp;
  if($_=~/\>/){
    $toggle = 0;
    my($pham_to_check)=($_=~/\>.+\_UID\_\d+\s+.*?\s+(\d+?)\s+\w+\s+\d+$/);
    if(!defined $pham_to_check){
      ($pham_to_check)=($_=~/\>.+\_UID\_\d+\s+.*?\s+(\d+?)\s+/);
      if(!defined $pham_to_check){
        die print "Could not find pham in $_!\n";
      }
    }
    foreach my $check (@phams){
      if($pham_to_check eq $check){
        $toggle=1;
        print "Found it!\nFrom $_\n\n";
        print OUT "$_\n";
        last;
      }
      else{
        next;
      }
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