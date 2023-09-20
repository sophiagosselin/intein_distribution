#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use Data::Dumper;

#23_09_18
#Sophia Gosselin
#takes a directory containing multiple fasta files w/ intein amino acid sequences
#these files must be the only files in said directory, but can be in a nested format.

#outputs the associated extein and complete acession amino acid sequences
#also outputs the associated nucleotide sequences for each associated type of sequence 
#(intein, extein, and complete acession)



#notes: need to make hash pairing everything from blast searches




#inputs/globals
#input 1 is location of fasta files
my $input_directory = $ARGV[0];

#input 2 is your protein database (ran with parse_seqids)
my $prot_database = $ARGV[1];

#input 3 is your nucleotide database (ran with parse_seqids)
my $nucl_database = $ARGV[2];


my(@intein_prot_fastas) = SETUP();
MAIN();

sub SETUP{
    #make directories, and get files needed for downstream processes
    my(@intein_fasta)=GET_INPUTS_RECURSIVE($input_directory);

    return(@intein_fasta);

}

sub MAIN{

    foreach my $fastafile (@intein_prot_fastas){

        #get complete protein sequences from intein sequences        
        my($full_prot_seq_file,$hashref) = GET_FULL_PROTEIN_SEQ($fastafile);
        my(%full_intein_paired_memory) = %{$hashref};

        #get extein sequence from complete sequence using intein sequences.
        my $extein_prot_seq_file = GET_EXTEIN_SEQ($full_prot_seq_file,$fastafile,\%full_intein_paired_memory);
        
        #GET_NUCLEOTIDES($fastafile);

    }
}

sub GET_INPUTS_RECURSIVE{
  #recrusively gets all inputs from a directory (even if in recursive directories)
  #only files with fasta faa or fna will be captured.
  my $directory = shift;
  my @files;

  if($directory!~/\/$/){
    $directory.="\/";
  }

  my @directory_files = glob "$directory*";
  
  #check if directories or not
  foreach my $input (@directory_files){
    if(-f "$input"){
      if("$input"=~/.*\.[fasta|faa|fna]/){
        push(@files,"$input");
      }
      else{
        next;
      }
    }
    else{
      my(@recursive_files)=GET_INPUTS_RECURSIVE("$input");
      push(@files,@recursive_files);
    }
  }
  return(@files);
}

sub GET_EXTEIN_SEQ{
    #takes intein and full sequence files as well as the memory hash that pairs ascessions
    #returns a file containing only the extein sequence of the full sequence
    my $full_seq_fasta = shift;
    my $intein_seq_fasta = shift;
    my %paried_memory = %{my $hashref = shift};
    my %full_seq_fastadata = READIN_FASTA($full_seq_fasta);
    my %intein_fastadata = READIN_FASTA($intein_seq_fasta);


    my($file_handle)=($full_seq_fasta=~/.*?\/(.*?)\.faa/);
    mkdir("extein_prot_seq_clusters");
    open(my $extein_out, "+> extein_prot_seq_clusters\/$file_handle\_extein.faa");
    foreach my $intein_asc (keys %intein_fastadata){
        #for each intein sequence, find the associated full sequence ascession
        my $intein_sequence = $intein_fastadata{$intein_asc};
        my $full_asc = $paired_memory{$intein_asc}{"full_seq_original_asc"};
        my $full_sequence = $full_seq_fastadata{$full_asc};

        #splice the intein out of the full sequence
        my($cterminal)=($full_sequence=~/^(.*?)$intein_sequence.*$/);
        my($nterminal)=($full_sequence=~/^.*?$intein_sequence(.*)$/);

        #put the pieces back together
        my $extein_sequence = $cterminal.$nterminal;
        
        #check for errors
        if(!defined $extein_sequence){
            die print "Could not extract the extein sequence for intein $intein_asc paired with the full sequence $full_asc\n
            For the purposes of debuggin, the intein sequence is $intein_sequence
            The full sequence is $full_sequence\n\n";
        }

        #to output
        print $extein_out "$full_asc\_extein\_only\n$extein_seq\n";
    }
    close $extein_out;

    return("extein_prot_seq_clusters\/$file_handle\_extein.faa");
}

sub GET_FULL_PROTEIN_SEQ{
    #takes fasta file of intein sequence, and file data hash as input.
    #retrieves the associated complete asc amino acid sequence
    my $fasta = shift;
    my(%fastadata)=READIN_FASTA($fasta);

    #blast file against database
    my($blast_file)=RUN_BLAST("blastp",$prot_database,$fasta,"\"6 qseqid sseqid pident\"");
    
    #parse blast file
    my($hashref2,$hashref3)=PARSE_BLAST($blast_file,\%fastadata);
    my %blast_matches = %{$hashref2};
    my %memory_hash=%{$hashref3};
    
    #check to see if you missed results
    my(@keys)=(keys %blast_matches);
    my $number_of_matches = @keys;
    my (@keys2) = (keys %fastadata);
    my $number_of_sequences = @keys2;
    if($number_of_matches != $number_of_sequences){
        die print "For file $fasta:\nWas expecting $number_of_sequences, but only found $number_of_matches\n";
    }

    #create range file and extract associated sequences
    mkdir("full_prot_seq_clusters");
    my($file_handle)=($fasta=~/.*\/(.*?)\./);
    EXTRACT_MATCHES(\%blast_matches,$prot_database,"full_prot_seq_clusters\/$file_handle.faa");

    #finaly pair up the extracted full seq ascs with their blast counterparts
    my $match_counter=0;
    my %full_seq_data = READIN_FASTA("full_prot_seq_clusters\/$file_handle.faa");
    foreach my $intein_asc (keys %memory_hash){
        #print "Searching for $intein_asc\n";
        foreach my $full_asc (keys %full_seq_data){
            my $check = $memory_hash{$intein_asc}{"full_seq_blast_asc"};
            #print "Comparing blast asc $check to full asc $full_asc\n";
            if($full_asc=~/$check/){
                #print "Found match!\n\n";
                $memory_hash{$intein_asc}{"full_seq_original_asc"}=$full_asc;
                $match_counter++;
                last;
            }
            else{
                #print "No match\n";
                next;
            }
        }
        if(!defined $memory_hash{$intein_asc}{"full_seq_original_asc"}){
            die print "\n\nCould not find match to $intein_asc in full seq ascs.\n\n";
        }
    }

    #check to see if all matches were paired properly
    if($match_counter != $number_of_matches){
        die print "\n\nMissing matches for full seq asc!\nWas expecting $number_of_matches, but only found $match_counter.\nCheck the logs above.\n\n";
    }

    return("full_prot_seq_clusters\/$file_handle.faa",\%memory_hash);
}

sub EXTRACT_MATCHES{
    #takes a hash of matches
    #extracts matches, print them to a range file and extract sequences
    #returns hash of sequences and ascessions
    my (%seqs_for_range)= %{my $hashref = shift};
    my $database = shift;
    my $output_name = shift;

    open(my $range, "+> range.txt");
    foreach my $asc (keys %seqs_for_range){
        my @split_data = split(/\t/,$seqs_for_range{$asc});
        print $range "$split_data[1]\n";
    }
    close $range;

    #extract matches
    system("blastdbcmd -db $database -entry_batch range.txt -outfmt \"%f\" > $output_name");
}

sub RUN_BLAST{
    #takes file to search, and whether the search is a tblastn or a blastp
    my ($blastversion) = shift;
    my ($database) = shift;
    my ($query_fasta) = shift;
    my ($out_version) = shift;
    system("$blastversion -query $query_fasta -db $database -out blast6.txt -evalue 1e-50 -outfmt $out_version");

    return("blast6.txt");
}

sub PARSE_BLAST{
    #blast file outfmt 6 as input
    #filters results to find query sequences in database
    #returns hash of matches (1 per query)
    my $blast_to_parse = shift;
    my %intein_data = %{my $hashref = shift};

    my %matches;
    my %memory;
    open(my $blast6, "< $blast_to_parse");
    while(<$blast6>){
        chomp;
        my @split = split(/\t/,$_);

        #percent ID cutoff
        if($split[2] ne "100.000"){
            #print "Percent ID of match only $split[4] not 100.\n";
            next;
        }

        #check if match has been previously found to the given query
        if(!defined $matches{$split[0]}){
            my($phagename)=GET_PHAGE_NAME($split[0]);
            if(!defined $phagename){
                die print "Failed to get phage name for $split[0]\n\n";
            }
            #check for matching name
            if($split[1]=~/$phagename/){
                $matches{$split[0]}=$_;
                my($paired_intein_match)=FIND_ASSOCIATED_FULL_ASC($split[0],\%intein_data);
                $memory{$paired_intein_match}{"intein_blast_asc"}=$split[0];
                $memory{$paired_intein_match}{"full_seq_blast_asc"}=$split[1];
                #print "$phagename matches $split[1]!\n";
            }
            
            #skip non match
            else{
                next;
            }

        }

        #skip match to already found query
        else{
            next;
        }
    
    }
    close $blast6;

    return(\%matches,\%memory);
}

sub FIND_ASSOCIATED_FULL_ASC{
    my $asc_to_find = shift;
    my %ascs = %{my $hashref = shift};
    my $found_it;
    
    foreach my $asc_in_file (keys %ascs){
        if($asc_in_file=~/$asc_to_find/){
            $found_it=$asc_in_file;
        }
        else{
            next;
        }
    }

    return($found_it);
}

sub GET_PHAGE_NAME{
    #takes blast query line as input
    #returns the name of the phage or dies

    my $text = shift;
    my $name = "";

    #try to match the phage name using one of the following patterns
    #match test 1
    ($name)=($text=~/.*?\_(.*?)\_.*/);
    if($name eq ""){
        #match test 2
        ($name)=($text=~/something/);
        if($name eq ""){
            #match test 3
            ($name)=($text=~/something else/);
        }

        else{
            die print "Could not extract phage name from $text\n";
        }
    }

    else{}

    return($name);
}

sub READIN_FASTA{
    #takes fasta as input
    #returns hash using ascessions for keys and sequences for values 
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