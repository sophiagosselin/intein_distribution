#!/usr/bin/perl -w
use strict;
use warnings;
use Math::Trig;
use Math::Trig qw(great_circle_distance deg2rad);

#TAKES directory of sequences files as input (make sure to denote the handle)
#TAKES table containing lattitude and longitude values for all sequences in the fomer directory
#RETURNS an associated file containing pairwise geographical distances for each sequence file

#USER INPUTS
my $seq_dir = $ARGV[0];
my $seq_handle = $ARGV[1];
my $lat_long_file = $ARGV[2];

MAIN();

sub MAIN{
    #readin lat long file as nested hash
    my(%lat_long_data)=READIN_TABLE($lat_long_file,1);
    
    #get list of sequence files
    my(@seq_files_in_dir)=glob("$seq_dir/*.$seq_handle");

    #for each array calculate pairwise distances
    foreach my $seq_file (@seq_files_in_dir){
        #readin file
        my(@acessions)=READIN_FASTA($seq_file);
        my %pairwise_distances;

        #open file for later use
        open(my $local_meta, "+> $seq_file.meta");

        #pairwise seq comparisons
        foreach my $seq1 (@acessions){
            foreach my $seq2 (@acessions){
                
                #if same 2 acessions - set distance to 0
                if($seq1 eq $seq2){
                    $pairwise_distances{$seq1}{$seq2}=0;
                }

                #otherwise calculated distance in kms
                else{
                    my($distance) = GET_KM_DISTANCE("$seq1\t$seq2",$lat_long_data{$seq1}{"long"},$lat_long_data{$seq1}{"lat"},$lat_long_data{$seq2}{"long"},$lat_long_data{$seq2}{"lat"});
                    $pairwise_distances{$seq1}{$seq2}=$distance;
                }
            }
            
            #also print lat and long for these specific seqs to subtable for easier downstream use
            #comment out if not needed
            print $local_meta "$seq1\t$lat_long_data{$seq1}{\"lat\"}\t$lat_long_data{$seq1}{\"long\"}\n";

        }
        close $local_meta;

        #convert pairwise hash to matrix in string form and print to file
        my($toprint)=MATRIX_FROM_HASH(\%pairwise_distances);
        open(my $out, "+> $seq_file.distances");
        print $out "$toprint";
        close $out;

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

    my $reporter = shift;
    my $loc1_long = shift;
    my $loc1_lat = shift;
    my $loc2_long = shift;
    my $loc2_lat = shift;

    my $loc1_long_rad = deg2rad($loc1_long);
    my $loc1_lat_rad = deg2rad($loc1_lat);
    my $loc2_long_rad = deg2rad($loc2_long);
    my $loc2_lat_rad = deg2rad($loc2_lat);

    my $distance = great_circle_distance($loc1_long_rad,$loc1_lat_rad,$loc2_long_rad,$loc2_lat_rad,6378.137);
    #6372.8 is radius of earth in km. SO results are measured in km's

    return($distance);
}


sub READIN_FASTA{
  #takes name of fasta formatted file as input
  #returns array of ascessions 
  my $filehandle = shift;
  my @ascs;
  open(my $in, "< $filehandle");
  while(<$in>){
    chomp;

    if($_=~/\>/){
      push(@ascs,$_);
    }
    else{
      next;
    }

  }
  close $in;
  return(@ascs);
}

sub READIN_TABLE{
    #reads in a table file with or w/o header
    #column 1 is the entry name, #2 is the lattitude #3 is the longitude
    my $infile = shift;
    my $header_toggle = shift; #0 for no header, 1 for header
    my %output_hash;

    open(my $in, "< $infile");
    while(<$in>){
        chomp;

        #skips the first line if a header
        if($header_toggle == 1){
            $header_toggle = 0;
            next;
        }
        else{
            my(@data) = split(/\t/,$_);
            my($key)=$data[0];
            
            $output_hash{$key}{"lat"}=$data[1];
            $output_hash{$key}{"long"}=$data[2];
        }
    }
    close $in;
    return(%output_hash);
}