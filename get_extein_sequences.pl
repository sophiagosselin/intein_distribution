#!/usr/bin/perl
use strict;
use warnings;
#can  be commented out
use File::Copy;
use Data::Dumper;

#23_09_06
#Sophia Gosselin
#takes a directory of multiple fasta files containing intein sequences
#outputs the associated extein sequence.

#inputs/globals
#input 1 is location of fasta files
my $input_directory = $ARGV[0];
if($input_directory!~/\//){
  $input_directory.="\/";
}
#input 2 is your protein database (ran with parse_seqids)
my $database = $ARGV[1];

MAIN();

sub MAIN{

  mkdir("extein_prot_seq_clusters");
  mkdir("full_prot_seq_clusters");

  my @input_files = glob "$input_directory*";

  foreach my $fastafile (@input_files){
    #these lines can be commented out. For my use only.
    #my($new_dir_name)=($fastafile=~/.*\/(.*?)\..*/);
    #mkdir($new_dir_name);

    my($hashmem,$full_seq_file)=GET_FULL_SEQUENCE($fastafile);
    GET_EXTEIN_ONLY($fastafile,$full_seq_file,$hashmem);

    #same as above comment
    #my($no_dir)=($fastafile=~/.*?\/(.*?\.fasta)/);
    #copy($fastafile,"$new_dir_name\/$no_dir");
  }

}

sub GET_FULL_SEQUENCE{
  my $input_fasta = shift;

  #readin
  my(%paired_sequences)=READIN_FASTA($input_fasta);

  #blastp
  print "Running blastp for $input_fasta\n\n";
  system("blastp -query $input_fasta -db $database -out blast6.txt -evalue 1e-50 -outfmt \"6 qseqid sseqid pident\"");

  #parse BLAST output
  my $strand;
  my $counter=0;
  my %no_match;
  my %memory;
  open(my $range, "+> range.txt");
  open(my $blast6, "< blast6.txt");
  while(<$blast6>){
    chomp;
    my @split = split(/\t/,$_);

    #percent ID cutoff
    if($split[2] ne "100.000"){
      #print "Percent ID of match only $split[4] not 100.\n";
      next;
    }

    #you will need to change this line based on your database, and what the matches are reported as.
    #this pattern needs to match something unique to each query sequence such that it can be backtraced.
    my ($subject_name) = ($split[1]=~/(.*?\_\d+?)\_.*/);
    if(!defined $subject_name){
      ($subject_name)=($split[1]=~/(.*?\_.+?)\_.*/);
      if(!defined $subject_name){
        print "Could not extract name of subject sequence. Dying now.\n\nProblem line: $_\n";
        die;
      }
    }
    else{}
      
    if($split[0]=~/$subject_name/){
      #reporting for testing
      #print "Found match to $split[0].\nFrom match $_\n";

      print $range "$split[1]\n";
      if(!defined $memory{"$split[0]"}){
        $memory{"$split[0]"}="$split[1]";
      }
      elsif($memory{"$split[0]"} eq "0"){
        $memory{"$split[0]"}="$split[1]";
      }
      else{}
    }

    else{
        $no_match{"$split[0]"}{$counter}=$_;
        $counter++;
        if(!defined $memory{"$split[0]"}){
          $memory{"$split[0]"}="0";
          next;
        }
        elsif($memory{"$split[0]"} ne "0"){
          #print "\n\nNo match found to: $_\n\n";
          next;
        }
        else{
          $memory{"$split[0]"}="0";
          next;
        }
    }
  }
  close $blast6;

  foreach my $check (keys %memory){
    if($memory{$check} ne "0"){
      next;
    }
    else{
      foreach my $potential (keys %no_match){
        foreach my $inc (keys %{$no_match{$potential}}){
          my @data = split(/\t/,$no_match{$potential}{$inc});
          my ($phagename)=($data[1]=~/(.*?)\_\d+?\_.*/);
          if($potential=~/$phagename/){
            if($data[2] eq "100.000"){
              print $range "$data[1]\n";
              $memory{$check}=$data[1];
              last;
            }
            else{
              next;
            }
          }
          else{
            next;
          }
        }
      }
    }
  }
  close $range;
  
  #get full seqs
  print "Extracting full protein sequences\n\n";
  my($file_handle)=($input_fasta=~/.*?\/(.*?)\.fasta/);
  mkdir("full_prot_seq_clusters\/$file_handle");
  system("blastdbcmd -db $database -entry_batch range.txt -outfmt \"%f\" > full_prot_seq_clusters\/$file_handle\/$file_handle.faa");

  return(\%memory,"full_prot_seq_clusters\/$file_handle\/$file_handle.faa");
}

sub GET_EXTEIN_ONLY{
  my $intein_fasta = shift;
  my $full_seq_file = shift;
  my(%paired_mem)=%{my $hashmem = shift};

  my(%intein_content)=READIN_FASTA($intein_fasta);
  my(%full_seq_content)=READIN_FASTA($full_seq_file);

  my($file_handle)=($intein_fasta=~/.*?\/(.*?)\.fasta/);
  mkdir("extein_prot_seq_clusters\/$file_handle");
  open(my $extein, "+> extein_prot_seq_clusters\/$file_handle\/$file_handle.faa");

  foreach my $key (keys %full_seq_content){
    
    #get assoicated intein asc
    my($partial_key)=($key=~/^\>(.*?\_\d+?)\_.*/);
    if(!defined $partial_key){
      ($partial_key)=($key=~/(.*?\_.+?)\_.*/);
      if(!defined $partial_key){
        print "MOTHERFUCKER $key\n";
        die;
      }
    }  
    #print "Partial Key $partial_key\n";
    my($complete_int_key) = "";
    foreach my $int_key (keys %intein_content){
      #print "Int Key $int_key\n\n";
      if($int_key=~/$partial_key/){
        $complete_int_key = $int_key;
      }
      else{
        next;
      }
    }
    if($complete_int_key eq ""){
      foreach my $tempkey (keys %paired_mem){
        if($paired_mem{$tempkey}=~/$partial_key/){
          $complete_int_key=">$tempkey";
        }
        else{
          next;
        }
      }
    }

    if($complete_int_key eq ""){
      print "I HATE THIS SHIT $KEY\n";
      die;
    }
    
    my $intein_seq = $intein_content{$complete_int_key};
    my $full_seq = $full_seq_content{$key};

    my($cterminal)=($full_seq=~/^(.*?)$intein_seq.*$/);
    my($nterminal)=($full_seq=~/^.*?$intein_seq(.*)$/);

    my $extein_seq = $cterminal.$nterminal;
    
    print $extein "$key\n$extein_seq\n";
  }
  close $extein;
}

sub READIN_FASTA{
    #also standardizes any file it reads in
    my $infile = shift;
    my $accession="";
    my %sequences;
    
    open(IN, "< $infile");
    while(<IN>){
        chomp;
        $_=~s/[\/\-]/\_/g;
        if($_=~/\>/){
            $accession=$_;
            $sequences{$accession}="";
        }
        else{
            $sequences{$accession}.=$_;
        }
    }
    close IN;
    return(%sequences);
}

sub STANDARDIZE {
	#removes most unique characters from annotation lines
	#makes later searches and moving of files much easier.
	my $fastafile = shift;
	open(IN, "< $fastafile");
	open(OUT, "+> temp.fasta");
	while(<IN>){
		if($_=~/\>/){
			$_=~s/[\ \[\]\(\)\:\;\/\.\-]/\_/g;
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
