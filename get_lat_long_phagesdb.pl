#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use threads::shared;
use Thread::Queue;
use LWP::Simple;

#GLOBALS
#CHANGE PARAMETERS HERE IF DESIRED
my $input_file = $ARGV[0]; #input file first!
my $max_threads = $ARGV[1]; #check your number of cpus!

#open global files and thread queue
open(my $error_log, "+> error.log");
open(my $meta, "+> all_metadata.table");
print $meta "phage_name\tlattitude_choordinate\tlongitude_choordinate\n";
my $queue = Thread::Queue->new();

MAIN();
close $error_log;
close $meta;

sub MAIN{
  my (%input_sequences)=READIN_FASTA($input_file);

  my @input_ascessions=( keys %input_sequences );

  my (%long_lat)=GET_METADATA_PARENT(@input_ascessions);

  foreach my $phage_key (keys %long_lat){
    my $lat_toprint = $long_lat{$phage_key}{"lat"};
    my $long_toprint = $long_lat{$phage_key}{"long"};
    print $meta "$phage_key\t$lat_toprint\t$long_toprint\n";
  }
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
    my($phagename)=($asc=~/\>(.*?)\_.*/);
    #print "Captured the phage name $phagename!\n";
    push(@phagenames, $phagename);
    $paired_keys{$asc}=$phagename;
    $count_reporter++;
  }
  print "There are $count_reporter sequences to work through!\nStarting threads now\n";

  my($array_of_hashes_ref) = THREADED_TASK("GET_METADATA_CHILD",\@phagenames);
  my(@array_of_hashes) = @{$array_of_hashes_ref};
	my(%unparsed_metadata) = MERGE_NESTED_HASHES(@array_of_hashes);
  print "All threads finished\n";

  my %parsed_metadata;
  foreach my $key (keys %paired_keys){
    my $phagekey = $paired_keys{$key};
    $parsed_metadata{$key}{"lat"}=$unparsed_metadata{$phagekey}{"lat"};
    $parsed_metadata{$key}{"long"}=$unparsed_metadata{$phagekey}{"long"};
  }

  return(%parsed_metadata);
}
sub GET_METADATA_CHILD {
  #retrieves metadata for 1 phage
  #returns latitude and longitude values
  my %geo_data;
  while (my $phage = $queue->dequeue()) {
    my $url = "https://phagesdb.org/api/phages/$phage/";
    my $response = get($url);
    my($lat_choord,$long_choord);

    if (!$response) {
      print $error_log "Failed to retrieve metadata for $phage\n\n";
      $lat_choord="NA";
      $long_choord="NA";
    }

    else{
      my ($lat) = ($response =~ /\"found_GPS_lat\"\:.*?\"(.*?)\"/) or die print $error_log "Could not find lattitude in the following metadata!\n\n$response\n\n";;
      my ($long) = ($response =~ /\"found_GPS_long\"\:.*?\"(.*?)\"/) or die print $error_log "Could not find longitude in the following metadata!\n\n$response\n\n";
      
      #check if some dumbass inputted their data into incorrect fields.
      if($lat =~ /[WE]/ | $long =~ /[NS]/){
        $long_choord = STANDARDIZE_CHOORDS($lat,0);
        $lat_choord = STANDARDIZE_CHOORDS($long,1);
      }
      #normal pathway here
      else{
        $lat_choord = STANDARDIZE_CHOORDS($lat,1);
        $long_choord = STANDARDIZE_CHOORDS($long,0);
      }
    }

    #save choords and return
    $geo_data{$phage}{"lat"} = "$lat_choord";
    $geo_data{$phage}{"long"} = "$long_choord";
  }
  return \%geo_data;
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
  return(\@return_values);
}

sub MERGE_NESTED_HASHES{
	my @hashes_to_merge = @_;
	my %new_hash;
	foreach my $hashref (@hashes_to_merge){
		my %temp_hash = %{$hashref};
		%new_hash = (%new_hash,%temp_hash);
	}
	my $test = (\%new_hash);
	return(%new_hash);
}

sub READIN_FASTA{
  #takes name of fasta formatted file as input
  #returns hash of sequences with ascessions as keys
  my $filehandle = shift;
  my %seq_data;
  my $asc_holder="";
  open(IN, "< $filehandle");
  while(<IN>){
    chomp;
    if($_=~/\>/){
      $asc_holder=$_;
      $seq_data{$asc_holder}="";
    }
    else{
      $seq_data{$asc_holder}.=$_;
    }
  }
  close IN;
  return(%seq_data);
}

sub STANDARDIZE_FASTA {
	#removes most unique characters from annotation lines
	#makes later searches and moving of files much easier.
  #also tacks a uniqueid on the end
	my $fastafile = shift;
  my $counter = 0;
	open(IN, "< $fastafile");
	open(OUT, "+> temp.fasta");
	while(<IN>){
    chomp;
		if($_=~/\>/){
			$_=~s/[\ \[\]\(\)\:\;\/\.\-\~\`\!\@\#\$\%\^\&\*\=\+\{\}\?]/\_/g;
			print OUT "$_\_$counter\n";
      $counter++;
		}
		else{
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;
	unlink $fastafile;
	rename "temp.fasta", $fastafile;
}

sub STANDARDIZE_CHOORDS{
  my $choord = shift;
  my $lat_or_long = shift; #1 for lat, 0 for long
  my $choord_standardized = "";

  if($lat_or_long == 1){
    #it's lattitude! Do stuff here
    if($choord =~ /[NS]/){
      my($lat,$direction) = ($choord =~/(.*?)([NS])/);
      if($direction eq "N"){
        $choord_standardized = $lat;
      }
      elsif($direction eq "S"){
        $choord_standardized = ($lat*-1);
      }
      else{
        die "Error. Check $choord. Found $lat and $direction\n";
      }
    }
    else{
      $choord_standardized = $choord; 
    }
  }
  elsif($lat_or_long == 0){
    #it's longitude!
    if($choord =~ /[WE]/){
      my($long,$direction) = ($choord =~/(.*?)([WE])/);
      if($direction eq "E"){
        $choord_standardized = $long;
      }
      elsif($direction eq "W"){
        $choord_standardized = ($long*-1);
      }
      else{
        die "Error. Check $choord. Found $long and $direction\n";
      }
    }
    else{
      $choord_standardized = $choord; 
    }
  }

  chomp($choord_standardized);
  return($choord_standardized);
}