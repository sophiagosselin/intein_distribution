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

#FUTURE DIRECTIONS
#find any clusters where there is a distance value above X threshold
#make phylogenies of those clusters.

#compare clustering at 98 98.5 99 and 99.5% seqid

#send peter info on how much money you were promised over the summer from any sources.

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
my $cluster_id = $ARGV[1]; #suggest .985 or .99
my $max_threads = $ARGV[2]; #check your number of cpus!
my $blast_database = $ARGV[3]; #path to blastdb that has seqs of interest. must be made using parse_seqids

MAIN();
close $error_log;
close $meta;

sub MAIN{
  my (%input_sequences)=READIN_FASTA($input_file);

  my @input_ascessions=( keys %input_sequences );

  my (%long_lat)=GET_METADATA_PARENT(@input_ascessions);

  my (@cluster_fastas)=USEARCH($input_file,$cluster_id);

  open(my $AVG, "> average_distances.table");
  print $AVG "Cluster_Number\tNumber_of_Members\tAverage_Distance_km\tUpper_Distance_Outlier\tExtein_Sequence_Similarity\tLower_Similarity_Outlier\n";

  #note that we skip any clusters with less than 2 taxa
  foreach my $cluster (@cluster_fastas){
    my(%cluster_sequences)=READIN_FASTA($cluster);
    my @cluster_ascessions = ( keys %cluster_sequences );
    my $cluster_size = scalar @cluster_ascessions;
    my($clust_num)=($cluster=~/clusters\/(.*)/);
    if(scalar(@cluster_ascessions)>=3){
      my %cluster_distances;
      foreach my $asc1 (@cluster_ascessions){
        foreach my $asc2 (@cluster_ascessions){
          my $distance = GET_KM_DISTANCE($long_lat{$asc1}{"long"},$long_lat{$asc1}{"lat"},$long_lat{$asc2}{"long"},$long_lat{$asc2}{"lat"},$asc1,$asc2);
          $cluster_distances{$asc1}{$asc2}=$distance;
        }
      }
      #print out matrix for later supplemental data
      my $toprint = MATRIX_FROM_HASH(\%cluster_distances);
      open(OUT, "+> $cluster.dist");
      print OUT $toprint;
      close OUT;

      #now get data for the amount of sequence divergence of the associated exteins
      my ($similarity_data,$similarity_outlier_up,$similarity_outlier_low)=EXTEIN_DIVERGENCE($clust_num,\%cluster_sequences);
      #calculate average distance between all members of the cluster
      #note that any cell with NA is excluded.
      my $sum=0;
      my $count=0;
      my $upper_outlier="NA";
      my $lower_outlier="NA";
      foreach my $row_asc (keys %cluster_distances){
        foreach my $col_asc (keys %{$cluster_distances{$row_asc}}){
          if($cluster_distances{$row_asc}{$col_asc} eq "NA"){
            next;
          }
          else{
            $sum+=$cluster_distances{$row_asc}{$col_asc};
            $count++;

            if($upper_outlier eq "NA"){
              $upper_outlier=$cluster_distances{$row_asc}{$col_asc};
            }
            elsif($upper_outlier<$cluster_distances{$row_asc}{$col_asc}){
              $upper_outlier=$cluster_distances{$row_asc}{$col_asc};
            }
            else{}

            if($lower_outlier eq "NA"){
              $lower_outlier=$cluster_distances{$row_asc}{$col_asc};
            }
            elsif($lower_outlier>$cluster_distances{$row_asc}{$col_asc}){
              $lower_outlier=$cluster_distances{$row_asc}{$col_asc};
            }
            else{}
          }
        }
      }
      if($count == 0){
        print $AVG "$clust_num\t$cluster_size\tINSUFFICIENT_GEO_DATA\t$upper_outlier\t$similarity_data\t$similarity_outlier_low\n";
      }
      else{
        my $average = $sum/$count;
        print $AVG "$clust_num\t$cluster_size\t$average\t$upper_outlier\t$similarity_data\t$similarity_outlier_low\n";
      }
    }
    else{
      #TFTC too few to count
      print $AVG "$clust_num\t$cluster_size\tNA\tNA\tNOT_ENOUGH_MEMBERS_FOR_CALCULATIONS\tNA\n";
    }
  }
}

sub EXTEIN_DIVERGENCE{
  #takes hash of sequences and associated ascessions
  #runs blast searches using said sequences against a blastdb protein database
  #returns the full sequence from the best match
  my $cluster_number=shift;
  my $hashref=shift;
  my %intein_sequences = %{$hashref};
  mkdir("blast_searches");

  open(EXT, "+> blast_searches/$cluster_number\_extein.fasta");
  foreach my $accession (keys %intein_sequences){

    my($filename,$phagename)=($accession=~/\>(.*?\_(.*?\_.*?)\_).*/);

    open(OUT, "+> blast_searches\/$filename\.fasta");
    print OUT "$accession\n";
    print OUT "$intein_sequences{$accession}\n";
    close OUT;

    system("blastp -db $blast_database -query blast_searches/$filename.fasta -evalue 1e-50 -outfmt 6 -out 'blast_searches/$filename.blast'");
    if(-s "blast_searches/$filename.blast"){
      FILTER_BLAST($phagename,"blast_searches/$filename.blast");
      system("blastdbcmd -db $blast_database -entry_batch blast_searches/filtered_results.txt -outfmt \"\%f\" > blast_searches/$filename\_extein.fasta");
      if(-z "blast_searches/$filename\_extein.fasta"){
        print $error_log "$accession could not be extracted. Testing for acession changes.\n";
        my($no_cds)=($phagename=~/(.*?)\_.*/);
        FILTER_BLAST($no_cds,"blast_searches/$filename.blast");
        system("blastdbcmd -db $blast_database -entry_batch blast_searches/filtered_results.txt -outfmt \"\%f\" > blast_searches/$filename\_extein.fasta");
        if(-z "blast_searches/$filename\_extein.fasta"){
          print $error_log "$accession could not be extracted even after ignoring cds. Excluding it from sequence similarity calculations.\n";
          next;
        }
        else{}
      }
      else{}
      my %extein_sequences = READIN_FASTA("blast_searches/$filename\_extein.fasta");
      foreach my $ext_asc (keys %extein_sequences){
        my $ext_seq_w_int = $extein_sequences{$ext_asc};
        my $int_seq = $intein_sequences{$accession};
        $ext_seq_w_int =~ s/$int_seq//g;
        print EXT "$ext_asc\n";
        print EXT "$ext_seq_w_int\n";
      }
    }
    else{
      print $error_log "$accession returned no matches from BLAST. Excluding it for the purposes of calculaing sequence similarity.\n";
    }
  }
  close EXT;

  STANDARDIZE_FASTA("blast_searches/$cluster_number\_extein.fasta");
  system("muscle -align blast_searches/$cluster_number\_extein.fasta -output blast_searches/$cluster_number\_extein.fasta.aligned");

  my($sequence_similarity,$seq_out_up,$seq_out_low) = SEQUENCE_SIMILARITY("blast_searches/$cluster_number\_extein.fasta.aligned");

  return($sequence_similarity,$seq_out_up,$seq_out_low);
}

sub SEQUENCE_SIMILARITY{
  #takes alignment as input
  #calculates and returns the average sequence similarity
  my $sequence_file = shift;
  my %sequence_data = READIN_FASTA($sequence_file);
  my @pairwise_seq_sim;

  my $sim_out_up = "NA";
  my $sim_out_low = "NA";
  foreach my $seq_asc1 (keys %sequence_data){
    my $sequence_length = length($sequence_data{$seq_asc1});
    #print "Length of $seq_asc1 is $sequence_length\n";
    foreach my $seq_asc2 (keys %sequence_data){
      my $similarity_counter = 0;
      my $length_to_average_over = $sequence_length;
      my @split_seq2 = split(//,$sequence_data{$seq_asc2});
      my @split_seq1 = split(//,$sequence_data{$seq_asc1});
      for (my $index=0; $index<$sequence_length; $index++){
        if($split_seq1[$index] eq "-" | $split_seq2[$index] eq "-"){
          $length_to_average_over--;
        }
        elsif($split_seq1[$index] eq $split_seq2[$index]){
          $similarity_counter++;
        }
        else{
          #Nothing!
        }
      }
      my $similarity = ($similarity_counter/$length_to_average_over);

      if($sim_out_low eq "NA"){
        $sim_out_low=$similarity ;
      }
      elsif($sim_out_low>$similarity ){
        $sim_out_low=$similarity ;
      }
      else{}

      if($sim_out_up eq "NA"){
        $sim_out_up=$similarity ;
      }
      elsif($sim_out_up<$similarity ){
        $sim_out_up=$similarity ;
      }
      else{}
      #print "Similarity is $similarity.\nTotal Length of Alignment sans uninformative sites is $length_to_average_over\nSimilar Sites count is $similarity_counter\n\n";
      push(@pairwise_seq_sim,$similarity);
    }
  }

  #calculate average percent similarity
  my $arr_size = @pairwise_seq_sim;
  my $running_average = 0;
  foreach my $percent (@pairwise_seq_sim){
    $running_average+=$percent;
  }
  my $avg_sim = $running_average/$arr_size;
  return($avg_sim,$sim_out_up,$sim_out_low);
}

sub FILTER_BLAST{
	#takes a blast output file (format 6) as input.
	#Retrieves the best hits for each match so they can be extracted.
  my $phage_of_interest=shift;
  my $blast_results = shift;
  my $toggle = 0;
	my @entries_for_blastcmd;
	open(BLAST, "< $blast_results") or die VERBOSEPRINT(0, "Check your BLAST software and databases for issues. No BLAST output from search was found.\n");
	open(OUT, "> blast_searches/filtered_results.txt");
	while(<BLAST>){
		chomp;
    my @output_columns = split(/\t/, $_);
    if($output_columns[1]=~/$phage_of_interest/){
      if($toggle==0){
        print OUT "$output_columns[1]\n";
        $toggle=1;
      }
      elsif($toggle==1){
        #print "Found $phage_of_interest twice... in $blast_results.\n";
      }
      else{
        next;
      }
    }
    else{
      next;
    }
	}
	close BLAST;
	close OUT;
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
  my $asc_test1 = shift;
  my $asc_test2 = shift;

  my @test_arr=($loc1_lat,$loc1_long,$loc2_lat,$loc2_long);
  foreach my $test (@test_arr){
    if($test eq "NA"){
      return("NA");
    }
    elsif($test eq ""){
      print "Unititialized values from $asc_test1 and $asc_test2\n";
    }
    else{
      next;
    }
  }

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
    print $meta "$phage has been parsed.\n";
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
      if(!$lat){
        print $error_log "Could not find lattitude in the following metadata!\n\n$response\n\n";
        $lat_choord="NA";
      }
      elsif($lat=~/.*([NS])/){
        my $direction = $1;
        #check this later.
        ($lat_choord)=($lat=~/(.*?)\ .*/);
        if(!defined $lat_choord){
          $lat_choord="NA";
        }
        elsif($direction eq "N"){
          $lat_choord=(90-$lat_choord);
        }
        elsif($direction eq "S"){
          $lat_choord=(abs($lat_choord)+90);
        }
        else{
          print $error_log "Program claimed lattitude: $lat, and claimed the direction was $direction.\nThis should contain a N or S, but neither was captured.\nAssociated metadata is:$response\n\n";
          $lat_choord="NA";
        }
      }
      else{
        if($lat=~/[WE]/){
          print $error_log "Program captured lattitude: $lat but could this is actually a longitude value.\nAssociated metadata is:$response\n\n";
          $lat_choord="NA";
        }
        elsif($lat=~/\-/){
          $lat_choord=(abs($lat)+90);
        }
        elsif($lat=~/\d+/){
          $lat_choord=(90-$lat);
          if(!defined $lat_choord){
            print $error_log "Program captured lattitude: $lat but could not capture the choord.\nAssociated metadata is:$response\n\n";
            $lat_choord="NA";
          }
        }
        else{
          print $error_log "Program captured lattitude: $lat but this is not a numeric value.\nAssociated metadata is:$response\n\n";
          $lat_choord="NA";
        }
      }

      my ($long) = ($response =~ /\"found_GPS_long\"\:.*?\"(.*?)\"/) or die print $error_log "Could not find longitude in the following metadata!\n\n$response\n\n";
      if(!$long){
        print $error_log "Could not find longitude in the following metadata!\n\n$response\n\n";
        $long_choord="NA";
      }
      elsif($long=~/([WE])/){
        my $direction = $1;
        ($long_choord)=($long=~/(.*?)\ .*/);
        if(!defined $long_choord){
          $long_choord="NA";
        }
        elsif($direction eq "E"){
          #nothing!
        }
        elsif($direction eq "W"){
          $long_choord=(abs($long_choord)+180);
        }
        else{
          print $error_log "Program claimed longitude: $long and claimed the direction was $direction.\nThis should contain a E or W, but neither was captured.\nAssociated metadata is:$response\n\n";
          $long_choord="NA";
        }
      }
      else{
        if($long=~/[NS]/){
          print $error_log "Program captured longitude: $long but could this is actually a lattitude value.\nAssociated metadata is:$response\n\n";
          $long_choord="NA";
        }
        elsif($long=~/\-/){
          $long_choord=(abs($long)+180);
        }
        elsif($long=~/\d+/){
          #the choordinate should be a choordinate
          $long_choord=$long;
          if($long_choord eq ""){
            print $error_log "Program captured lattitude: $long but could not capture the choord.\nAssociated metadata is:$response\n\n";
            $lat_choord="NA";
          }
        }
        else{
          print $error_log "Program captured longitude: $long but this is not a numeric value.\nAssociated metadata is:$response\n\n";
          $long_choord="NA";
        }
      }

      if($lat_choord eq ""){
        $lat_choord="NA";
      }
      else{}
      if($long_choord eq ""){
        $long_choord="NA";
      }
      else{}
    }
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
