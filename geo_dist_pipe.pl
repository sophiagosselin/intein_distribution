#!/usr/bin/perl -w
use strict;
use warnings;
use threads;
use threads::shared;
use Thread::Queue;
use LWP::Simple;
use Math::Trig;
use Math::Trig qw(great_circle_distance deg2rad);

#FUTURE DIRECTIONS
#find any clusters where there is a distance value above X threshold
#make phylogenies of those clusters.

#CLEANUP AND ERROR PREP
system("rm error.log");
system("rm *.temp");
system("rm -r clusters");
system("rm average_distances.table");
system("rm clustered.uc");
system("rm average_distances.table");
open(my $error_log, "+> error.log");
my $queue = Thread::Queue->new();

#GLOBALS
#CHANGE PARAMETERS HERE IF DESIRED
my $input_file = $ARGV[0]; #input file first!
my $cluster_id = $ARGV[1]; #suggest .985 or .99
my $max_threads = $ARGV[2]; #check your number of cpus!

MAIN();
close $error_log;

sub MAIN{
  my(%input_sequences)=READIN_FASTA($input_file);

  my @input_ascessions = ( keys %input_sequences );

  my(%long_lat)=GET_METADATA_PARENT(@input_ascessions);

  my (@cluster_fastas) = USEARCH($input_file,$cluster_id);

  open(my $AVG, "+> average_distances.table");
  foreach my $cluster (@cluster_fastas){
    my(%cluster_sequences)=READIN_FASTA($cluster);
    my @cluster_ascessions = ( keys %cluster_sequences );
    my %cluster_distances;
    foreach my $asc1 (@cluster_ascessions){
      foreach my $asc2 (@cluster_ascessions){
        my $distance  = GET_KM_DISTANCE($long_lat{$asc1}{"long"},$long_lat{$asc1}{"lat"},$long_lat{$asc2}{"long"},$long_lat{$asc2}{"lat"});
        $cluster_distances{$asc1}{$asc2}=$distance;
      }
    }
    #print out matrix for later supplemental data
    my $toprint = MATRIX_FROM_HASH(\%cluster_distances);
    open(OUT, "+> $cluster.dist");
    print OUT $toprint;
    close OUT;

    #calculate average distance between all members of the cluster
    my $sum=0;
    my $count=0;
    foreach my $row_asc (keys %cluster_distances){
      foreach my $col_asc (keys %{$cluster_distances{$row_asc}}){
        $sum+=$cluster_distances{$row_asc}{$col_asc};
        $count++;
      }
    }
    my $average = $sum/$count;
    my($clust_num)=($cluster=~/clusters\/(.*)/);
    print $AVG "$clust_num\t$average\n";
  }
}

sub MATRIX_FROM_HASH{
	#takes hash of data to convert into a 2D matrix, as well as an array of headers
	#returns the 2D matrix as a string
	my %matrix_data = %{my $hashref = shift};
	my $matrix_string ="";

	#create matrix header
  my @matrix_header = ( keys %matrix_data );
	foreach my $header (sort @matrix_header){
		$matrix_string.="$header\t";
	}

	#convert hash to matrix
	my $current_row = "";
	foreach my $row_entry (sort keys %matrix_data){
    foreach my $column_entry (sort keys %{$matrix_data{$row_entry}}){
      if($row_entry ne $current_row){
        $matrix_string.="\n$row_entry\t";
        $current_row = $row_entry;
      }
      my $distance_to_print = $matrix_data{$row_entry}{$column_entry};
      $matrix_string.="$distance_to_print\t";
		}
	}
	return($matrix_string);
}

sub GET_KM_DISTANCE{
  #takes long and lat for 2 points as inputs
  #returns distance in km between those two points

  my $loc1_long = shift;
  my $loc1_lat = shift;
  my $loc2_long = shift;
  my $loc2_lat = shift;

  my $loc1_long_rad = deg2rad($loc1_long);
  my $loc1_lat_rad = deg2rad($loc1_lat);
  my $loc2_long_rad = deg2rad($loc2_long);
  my $loc2_lat_rad = deg2rad($loc2_lat);

  my $distance = great_circle_distance($loc1_long_rad,$loc1_lat_rad, $loc2_long_rad,$loc2_lat_rad,6378.137);
  #6372.8 is radius of earth in km. SO results are measured in km's

  return($distance);
}

sub LAT_180_CONVERSION{
  # Notice the 90 - latitude: phi zero is at the North Pole.
  my $value = shift;
  my $converted_lat = abs(90-$value);
  return($converted_lat);
}

sub USEARCH{
  #takes fasta file and cluster % as input
  #runs usearch on file using provided parameters
  #returns all cluster files as an array

  my $fasta_file = shift;
  my $id = shift;
  mkdir("clusters");

  system("usearch -cluster_fast $fasta_file -sort length -id $id -uc clustered.uc -clusters clusters/cluster_");

  my @clusters = glob "clusters/*";
  foreach my $test (@clusters){
  }
  return(@clusters);
}

sub GET_METADATA_PARENT{
  #takes array of ascessions as inputs
  #returns metadata file associated with said ascessions
  my @ascs=@_;
  my @phagenames;
  my %paired_keys;
  my $count_reporter=0;

  foreach my $asc (@ascs){
    my($phagename)=($asc=~/\>.*?\_(.*?)\_.*/);
    push(@phagenames, $phagename);
    $paired_keys{$phagename}=$asc;
    $count_reporter++;
  }
  print "There are $count_reporter sequences to work through!\nStarting threads now\n";

  my($array_of_hashes_ref) = THREADED_TASK("GET_METADATA_CHILD",\@phagenames);
  my(@array_of_hashes) = @{$array_of_hashes_ref};
	my(%unparsed_metadata) = MERGE_NESTED_HASHES(@array_of_hashes);
  print "All threads finished\n";

  my %parsed_metadata;
  foreach my $key (keys %unparsed_metadata){
    my $newkey = $paired_keys{$key};
    $parsed_metadata{$newkey}{"lat"}=$unparsed_metadata{$key}{"lat"};
    $parsed_metadata{$newkey}{"long"}=$unparsed_metadata{$key}{"long"};
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

    if (!$response) {
      print $error_log "Failed to retrieve metadata for $phage\n";
    }

    else{
      my ($lat) = ($response =~ /"found_GPS_lat".*?\"(.*?)\".*[\,\}]/);
      if(!$lat){
        die print $error_log "Could not find lattitude in the following metadata!\n\n$response";
      }
      if($lat=~/([NS])/){
        my $dir = $1;
        ($lat)=($lat=~/(.*?)\ .*/);
        if($dir eq "N"){
          $lat=(90-$lat);
        }
        elsif($dir eq "S"){
          $lat=(abs($lat)+90);
        }
        else{
          die "Could not fix lattidue $lat";
        }
      }
      else{
        if($lat=~/\-/){
          $lat=(abs($lat)+90);
        }
        else{
          $lat=(90-$lat);
        }
      }
      my ($long) = ($response =~ /"found_GPS_long".*?\"(.*?)\".*[\,\}]/) or die print $error_log "Could not find longitude in the following metadata!\n\n$response";
      if($long=~/([WE])/){
        my $dir = $1;
        ($long)=($long=~/(.*?)\ .*/);
        if($dir eq "E"){
        }
        elsif($dir eq "W"){
          $long=(abs($long)+180);
        }
        else{
          die "Could not fix lattidue $long";
        }
      }
      else{
        if($long=~/\-/){
          $long=(abs($long)+180);
        }
        else{
        }
      }
      $geo_data{$phage}{"lat"} = $lat;
      $geo_data{$phage}{"long"} = $long;
    }
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
	my $fastafile = shift;
	open(IN, "< $fastafile");
	open(OUT, "+> temp.fasta");
	while(<IN>){
		if($_=~/\>/){
			$_=~s/[\ \[\]\(\)\:\;\/\.\-\~\`\!\@\#\$\%\^\&\*\=\+\{\}\?]/\_/g;
			print OUT $_;
		}
		else{
			print OUT $_;
		}
	}
	close IN;
	close OUT;
	unlink $fastafile;
	rename "temp.fasta", $fastafile;
}
