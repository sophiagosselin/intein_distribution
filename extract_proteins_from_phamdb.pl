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
open(INTERMIDIATE, "+> intermidiate.fasta");
while(<IN>){
  chomp;
  if($_=~/\>/){
    my($pham,$ascession) = ($_=~/\>\w+?\s+?\w+?\s+?(\w+?)\s+?(\w+?)\s+?\w+/);
    if(!$ascession){
      ($pham,$ascession) = ($_=~/\>\w+?\s+?\w+?\s+?(\w+?)\s+?(\w+)/);
    }
    if($pham eq "None"){
      next;
    }
    if($ascession eq $search_term){
      print INTERMIDIATE "$_\n";
    }
    else{
      next;
    }
  }
  else{
    next;
  }
}
close IN;
close INTERMIDIATE;

#takes intermidiate file and pulls the pham #
my @phams;
open(IN, "< intermidiate.fasta");
while(<IN>){
  chomp;
  my($pham) = ($_=~/\>.*?\ .*?\ (.*?)\ .*?\ .*/);
  push(@phams,$pham);
}
close IN;
unlink "intermidiate.fasta";

#get final extracted dataset
my $found_toggle = 0;
open(IN, "< $database");
open(OUT, "+> $search_term.fasta");
while(<IN>){
  chomp;
  if($_=~/\>/){
    $found_toggle = 0;
    my($pham) = ($_=~/\>\w+?\s+?\w+?\s+?(\w+?)\s+?\w+?\s+?\w+/);
        if(!$pham){
          ($pham) = ($_=~/\>\w+?\s+?\w+?\s+?(\w+?)\s+?\w+/);
        }
    if ( grep( /^$pham$/, @phams ) ) {
      print OUT "$_\n";
      $found_toggle = 1;
    }
    else{
      next;
    }
  }
  else{
    if($found_toggle == 1){
      print OUT "$_\n";
    }
    else{
      next;
    }
  }

}
close IN;
close OUT;
