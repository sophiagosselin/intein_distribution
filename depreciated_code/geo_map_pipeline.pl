#!/usr/bin/perl -w
use strict;
use warnings;

#USER DEFINED GLOBALS
my $input_seq = $ARGV[0];
my $thread_count = $ARGV[1];

#PIPE STARTS HERE
################################################################################################################3

#get lat and long data for all sequences in file
system("perl get_lat_long_phagesdb.pl $input_seq $thread_count");

#name of lat long table file
my $meta_file = "all_metadata.table";

#probably need to change the $input_seq her to the extein sequences of the inteins
system("usearch -cluster_fast $input_seq -sort length -id .975 -uc clustered.uc -clusters clusters/cluster_");

#add fasta extension to files
system("for f in clusters/*; do mv $f $f.fasta;done");

#get pairwise KM distance between seqs in clusters 
system("perl pairwise_geo_dist.pl clusters fasta $meta_file");

#take output and make a map figure for each cluster
