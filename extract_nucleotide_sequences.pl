#!/usr/bin/perl
use strict;
use warnings;
#can  be commented out
use File::Copy;

#23_09_06
#Sophia Gosselin
#takes a directory of multifasta files
#outputs all nucleotide sequences that 100% match said proteins
#this may lead to extra hits in cases of clonar genes


#inputs/globals
#input 1 is location of fasta files
my $input_directory = $ARGV[0];

#input 2 is your nucleotide database (ran with parse_seqids)
my $database = $ARGV[1];

#input 3 is name of output directory
my $output_directory = $ARGV[2];

my(@input_files) = GET_INPUTS_RECURSIVE($input_directory);

MAIN();

sub GET_INPUTS_RECURSIVE{
  #recrusively gets all inputs from a directory (even if in recursive directories)
  my $directory = shift;
  my @files;

  if($directory!~/\/$/){
    $directory.="\/";
  }

  my @directory_files = glob "$directory*";
  
  #check if directories or not
  foreach my $input (@directory_files){
    if(-f "$input"){
      push(@files,"$input");
    }
    else{
      my(@recursive_files)=GET_INPUTS_RECURSIVE("$input");
      push(@files,@recursive_files);
    }
  }
  return(@files);
}

sub MAIN{

  mkdir("$output_directory");

  foreach my $fastafile (@input_files){
    #these lines can be commented out. For my use only.
    #my($new_dir_name)=($fastafile=~/.*\/(.*?)\..*/);
    #mkdir($new_dir_name);

    GET_NUCLEOTIDES($fastafile);

    #same as above comment
    #my($no_dir)=($fastafile=~/.*?\/(.*?\.fasta)/);
    #copy($fastafile,"$new_dir_name\/$no_dir");
  }

}

sub GET_NUCLEOTIDES{
  my $input_fasta = shift;

  #readin
  my(%paired_sequences)=READIN_FASTA($input_fasta);

  #get seq length for coverage cuttoff
  my %cutoffs;
  foreach my $key (keys %paired_sequences){
    my $seq_length = length($paired_sequences{$key});
    my($key_for_blast)=($key=~/\>(.*?)\s.*/);
    $cutoffs{$key_for_blast}=($seq_length*3);
    print "Asc key from fasta is $key_for_blast\n";
  }

  #tblastn
  print "\n\nRunning tblastn for $input_fasta\n\n";
  system("tblastn -query $input_fasta -db $database -out blast6.txt -evalue 1e-50 -outfmt \"6 qseqid sseqid sstart send pident\"");

  #parse BLAST output
  my $strand;
  open(my $range, "+> range.txt");
  open(my $blast6, "< blast6.txt");
  while(<$blast6>){
    chomp;
    my @split = split(/\t/,$_);

    #percent ID cutoff
    if($split[4] ne "100.000"){
      #print "Percent ID of match only $split[4] not 100.\n";
      next;
    }
    else{}

    #you will need to change this line based on your database, and what the matches are reported as.
    #this pattern needs to match something unique to each query sequence such that it can be backtraced.
    my ($phage_name) = ($split[1]=~/(.*?)\_\d+$/);
    
    #check if you have captured the correct name.
    if(!defined $phage_name){
      #($subject_name)=($split[1]=~/.*?\|(.*?)\|.*/);
      if(!defined $phage_name){
        print "Could not extract name of subject sequence. Dying now.\n\nProblem line: $_\n";
        die;
      }
      else{}
      
    }
    
    #check if subject has correct phage name
    if($split[0]=~/$phage_name/){
      #reporting for testing
      #print "Found match to $split[0].\nFrom match $_\n";
      
      #coverage cutoff
      my $subject_length = abs($split[2]-$split[3])+1; #not sure why the plus one is needed here tbh
      my $backtrace_key = $split[0];
      
      #check if values are defined.
      if(!defined $subject_length){
        die print "Error at blast line $_\nSubject L $subject_length\nBacktrace $backtrace_key\nCutoff from Backtrace $cutoffs{$backtrace_key}\n\n";
      }
      if(!defined $cutoffs{$backtrace_key}){
        die print "Error at blast line $_\nSubject L $subject_length\nBacktrace $backtrace_key\nCutoff from Backtrace $cutoffs{$backtrace_key}\n\n";
      }
      if(!defined $backtrace_key){
        die print "Error at blast line $_\nSubject L $subject_length\nBacktrace $backtrace_key\nCutoff from Backtrace $cutoffs{$backtrace_key}\n\n";
      }

      #coverage cutoff
      if($subject_length != $cutoffs{$backtrace_key}){
        print "Length of $subject_length does not match $cutoffs{$backtrace_key}\n";
        next;
      }
      else{}

      #get values for range file
      if($subject_length >= 0){
        $strand = "minus";
        #print "Line is $_\n";
        print $range "$split[1]\ $split[3]\-$split[2]\ $strand\n";
      }
      else{
        #print "Line is $_\n";
        $strand = "minus";
        print $range "$split[1]\ $split[2]\-$split[3]\ $strand\n";
      }
    }

    else{
        #print "No match found.\n$_\n";
        next;
    }
  }
  close $blast6;
  close $range;
  
  #get nuc seqs
  print "Extracting nucleotide sequences\n\n";
  my($file_handle,$file_extension)=($input_fasta=~/.*\/(.*?)\.(fasta|faa)/);
  mkdir("$output_directory\/$file_handle");
  system("blastdbcmd -db $database -entry_batch range.txt -outfmt \"%f\" > $output_directory\/$file_handle\/$file_handle.fna");

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
