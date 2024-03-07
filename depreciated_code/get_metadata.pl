#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use threads::shared;
use Thread::Queue;
use LWP::Simple;
use Math::Trig;
use Math::Trig qw(great_circle_distance deg2rad);
use Data::Dumper;


#CLEANUP AND ERROR PREP
#system("rm error.log");
#system("rm meta.log");
#system("rm *.temp");
#system("rm -r clusters");
#system("rm -r blast_searches");
#system("rm average_distances.table");
#system("rm clustered.uc");
#system("rm average_distances.table");
open(my $error_log, "+> error.log");
open(my $meta, "+> meta.log");
my $queue = Thread::Queue->new();

#GLOBALS
#CHANGE PARAMETERS HERE IF DESIRED
my $input_file = $ARGV[0]; #input file first!
my $max_threads = $ARGV[1]; #check your number of cpus!

MAIN();

close $error_log;
close $meta;

sub MAIN{

    my (%input_sequences)=READIN_FASTA($input_file);

    my @input_ascessions=( keys %input_sequences );

    my(@toprint) = GET_METADATA_PARENT(@input_ascessions);

    open(my $out, "+> all_metadata.txt");
    foreach my $line (@toprint){
        print $out "$line";
    }
    close OUT;
}



sub GET_METADATA_PARENT{
  #takes array of ascessions as inputs
  #returns metadata file associated with said ascessions
  my @ascs=@_;
  my @phagenames;
  my %paired_keys;
  my $count_reporter=0;

  foreach my $asc (@ascs){
    #print "Acession is $asc!\n";
    my($phagename)=($asc=~/\>.*?\_(.*?)\_.*/);
    #print "Captured the phage name $phagename!\n";
    push(@phagenames, $phagename);
    $paired_keys{$asc}=$phagename;
    $count_reporter++;
  }
  print "There are $count_reporter sequences to work through!\nStarting threads now\n";

  my(@responses) = THREADED_TASK("GET_METADATA_CHILD",\@phagenames);
  print "All threads finished\n";

  return(@responses);
}
sub GET_METADATA_CHILD {
  #retrieves metadata for 1 phage
  #returns latitude and longitude values
  my %geo_data;
  while (my $phage = $queue->dequeue()) {
    print $meta "$phage has been parsed.\n";
    my $url = "https://phagesdb.org/api/phages/$phage/";
    my $response = get($url);
    my($lat_choord,$long_choord);

    if (!$response) {
        print $error_log "Failed to retrieve metadata for $phage\n\n";
        return("$phage\tNO DATA\n");
    }

    else{
        my($name)=($response=~/\"phage_name\"\:(.*?)[\,\}]/);
        my($cluster)=($response=~/\"cluster\"\:(.*?)[\,\}]/);
        my($subcluster)=($response=~/\"subcluster\"\:(.*?)[\,\}]/);
        my($city)=($response=~/\"found_city\"\:(.*?)[\,\}]/);
        my($state)=($response=~/\"found_state\"\:(.*?)[\,\}]/);
        my($country)=($response=~/\"found_country\"\:(.*?)[\,\}]/);
        my($lat)=($response=~/\"found_GPS_lat\"\:(.*?)[\,\}]/);
        my($long)=($response=~/\"found_GPS_long\"\:(.*?)[\,\}]/);
        return("$name\t$cluster\t$subcluster\t$city\t$state\t$country\t$lat\t$long\n");
    }
}

sub THREADED_TASK{
  my $subroutine = shift;
  my $array_ref = shift;
  my @thread_input = @{$array_ref};
  my @threads;
  my @return_values;

  foreach my $threadable_input (@thread_input){
    $queue->enqueue($threadable_input);
  }
  $queue->end();

  # Create a thread pool
  for (1..$max_threads) {
    push @threads, threads->create(\&$subroutine);
  }

  # Wait for all threads to finish
  foreach my $thread (@threads) {
    my $return_value = $thread->join();
    push(@return_values,$return_value);
  }
  return(@return_values);
}