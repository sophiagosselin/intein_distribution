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

#24_02_09
#Past Sophia wrote terrible code. Current Sophia hates you. AHHHHHH

#notes: need to make hash pairing everything from blast searches
#module load perl
#module load blast/2.11.0
#export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/allphams/
#export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/phagesdb/


#inputs/globals
#input 1 is location of fasta files
my $input_directory = $ARGV[0];

#input 2 is your protein database (ran with parse_seqids)
my $prot_database = $ARGV[1];

#input 3 is your nucleotide database (ran with parse_seqids)
my $nucl_database = $ARGV[2];

#START
#####################################################################

#READIN INPUT FILES
my(@intein_prot_fastas) = SETUP();
MAIN();

sub MAIN{

    mkdir("completed_inputs");

    foreach my $fastafile (@intein_prot_fastas){
        my %file_specific_memory_hash; #has following format: 
        #{intein aa ascs as keys from input file}
        ##{"full_aa_asc"} -> asc from file
        ##{"full_nuc_asc"} -> asc from file
        ##{"intein_nuc_asc"} -> asc from file

        #get complete protein sequences from intein sequences
        print "Getting full amino acid sequences for $fastafile\n";        
        my($full_prot_file,$hashref1,@queries_w_no_prot_match) = GET_FULL_PROTEIN_SEQ($fastafile);
        %file_specific_memory_hash = %{$hashref1};

        #get extein sequence from complete sequence using intein sequences.
        print "Getting extein amino acid sequences for $fastafile\n";  
        GET_EXTEIN_SEQ($full_prot_file,$fastafile,"extein_prot_seq","faa",0,\%file_specific_memory_hash,@queries_w_no_prot_match);
        
        print "Getting intein nucleotide sequences for $fastafile\n"; 
        my($intein_nuc_file,$array_ref1,$hashref2) = GET_NUCLEOTIDES($fastafile,1,"intein_nucl_seq",\%file_specific_memory_hash);
        my @queries_w_no_nuc_matches = @{$array_ref1};
        %file_specific_memory_hash = %{$hashref1};

        print "Getting full nucleotide sequences for $fastafile\n"; 
        my($full_nuc_file,$array_ref2,$hashref3) = GET_NUCLEOTIDES($full_prot_file,0,"full_nucl_seq",\%file_specific_memory_hash,@queries_w_no_nuc_matches);
        @queries_w_no_nuc_matches = @{$array_ref2};
        %file_specific_memory_hash = %{$hashref3};

        print "Getting extein nucleotide sequences for $fastafile\n"; 
        GET_EXTEIN_SEQ($full_nuc_file,$intein_nuc_file,"extein_nucl_seq","fna",1,\%file_specific_memory_hash,@queries_w_no_nuc_matches);

        #check if nucleotide files are missing content
        print "Completed work on $fastafile\n";
        my($copy_handle)=($fastafile=~/.*\/(.*)/);
        copy("$fastafile","completed_inputs\/$copy_handle");

        PRINT_NON_MATCHES(\@queries_w_no_prot_match,\@queries_w_no_nuc_matches,$fastafile);

        #for the sake of your sanity, this prints a key for which names are equivelant
        PRINT_KEY(\%file_specific_memory_hash);

    }
}

sub PRINT_KEY{
    my %key = %{ my $hashref1 = shift};

    #{intein aa ascs as keys from input file}
    ##{"full_aa_asc"} -> asc from file
    ##{"full_nuc_asc"} -> asc from file
    ##{"intein_nuc_asc"} -> asc from file


    open(my $keyfile, "+> accession_key.txt");
    print $keyfile "Intein_AA\tIntein_Nuc\tFull_Seq_AA\tFull_Seq_Nuc\tExtein_AA\tExtein_Nuc\n";

    foreach my $intein_aa_acc (keys %key){
        print $keyfile "$intein_aa_acc\t$key{$intein_aa_acc}{\"intein_nuc_asc\"}\t$key{$intein_aa_acc}{\"full_aa_asc\"}\t$key{$intein_aa_acc}{\"full_nuc_asc\"}\t$key{$intein_aa_acc}{\"full_aa_asc\"}\t$key{$intein_aa_acc}{\"full_nuc_asc\"}\n";
    }

    close $keyfile;
}

sub PRINT_NON_MATCHES{
    #prints non matches for given fasta file to new files
    my @prot_no_match = @{ my $arrayref = shift};
    my @nucl_no_match = @{ my $arrayref2 = shift};
    my $infasta = shift;
    my($file_handle)=($infasta=~/.*\/(.*?)\./);

    mkdir("no_matches");
    open(my $out1, "+> no_matches\/$file_handle\_no_matches_nucleotide\.txt");
    foreach my $no_match (@nucl_no_match){
        print $out1 "$no_match\n";
    }
    close $out1;

    open(my $out2, "+> no_matches\/$file_handle\_no_matches_aa\.txt");
    foreach my $no_match2 (@prot_no_match){
        print $out2 "$no_match2\n";
    }
    close $out2;
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

sub GET_NUCLEOTIDES{
    #takes protein fasta file as input 
    #returns the associated nucleotide sequence for each protein sequence in fasta file
    #saves acs to paired memory hash
    my $fasta_file = shift;
    my $mode = shift;
    my $outdir = shift;
    my %global_memory_copy = %{ my $hashref = shift };
    my @nomatch = @_;

    #blast file against database
    my($blast_file)=RUN_BLAST("tblastn -db_gencode 11 -seg no",$nucl_database,$fasta_file,"\"6 qseqid sseqid pident sstart send\"");

    #READIN query fasta
    my %fasta_data = READIN_FASTA($fasta_file);

    #parse BLAST file
    my($hashref2,$hashref3)=PARSE_BLAST($blast_file,\%fasta_data,1);
    my %blast_matches = %{$hashref2};
    my %local_memory_hash=%{$hashref3};

    #Extract matches to fna file
    mkdir("$outdir");
    move($blast_file,"$outdir\/$blast_file");
    my($file_handle)=($fasta_file=~/.*\/(.*?)\./);
    EXTRACT_MATCHES(\%blast_matches,$nucl_database,"$outdir\/$file_handle.fna",1);

    #append a unique ID to each line 
    APPEND_UID("$outdir\/$file_handle.fna");

    #check to see if you missed results. Push missed matches to a an array
    my %nucl_seq_fastadata = READIN_FASTA("$outdir\/$file_handle.fna");

    foreach my $prot_asc (keys %local_memory_hash){
        my $tracker = 0;

        foreach my $nucl_seq_asc (keys %nucl_seq_fastadata){

            #get values to confirm a match is paired correctly
            my $check_name = $local_memory_hash{$prot_asc}{"subject_blast_asc"};
            my $check_sstart = $local_memory_hash{$prot_asc}{"sstart"};
            my $check_sstop = $local_memory_hash{$prot_asc}{"send"};
            my @checkvals = ($check_name,$check_sstart,$check_sstop);
            my $toggle = 0;

            foreach my $check (@checkvals){
                #print "$nucl_seq_asc vs. $check\n";
                if($nucl_seq_asc=~/$check/){
                    next;
                }
                else{
                    $toggle = 1;
                }
            }

            #compare said asc to full seq asc 
            if($toggle == 0){
                #print "Nucl acs $nucl_seq_asc passed checkpoint 1!\n";
                #find the correct bin to place it in within the global memory hash
                foreach my $key1 (keys %global_memory_copy){
                    #check if this is an intein asc
                    if($mode == 1){
                        if($key1 eq $prot_asc){
                            $global_memory_copy{$prot_asc}{"intein_nuc_asc"}=$nucl_seq_asc;
                            #print "FOR DEBUGGING: The intein nuc asc for $prot_asc is $nucl_seq_asc\n";
                            #print "The BLAST results claim that $nucl_seq_asc equals $check_name\n\n";
                            $tracker = 1;
                            last;
                        }
                        else{
                            next;
                        }
                    }
                    else{
                        #check the full seq ascs
                        if($global_memory_copy{$key1}{"full_aa_asc"} eq $prot_asc){
                            $global_memory_copy{$key1}{"full_nuc_asc"}=$nucl_seq_asc;
                            #print "FOR DEBUGGING: The full nuc asc for $prot_asc is $nucl_seq_asc\n";
                            $tracker = 1;
                            last;
                        }
                        else{
                            next;
                        }
                    }
                }

                last;
            }

            else{
                next;
            }
        }

        #check if you could not find a match for the intein asc
        if($tracker == 0){
            print "\n\nCould not find match to $prot_asc in nuc seq ascs.\n";
            push(@nomatch,$prot_asc);
        }
    }

    return("$outdir\/$file_handle.fna",\@nomatch,\%global_memory_copy);
}

sub GET_EXTEIN_SEQ{
    #takes intein and full sequence files as well as the memory hash that pairs ascessions
    #returns a file containing only the extein sequence of the full sequence
    my $full_seq_fasta = shift;
    my $intein_seq_fasta = shift;
    my $dir_name = shift;
    my $extension = shift;
    my $flag = shift; #0 for aa, 1 for nucl
    my %paired_memory = %{my $hashref = shift};
    my @seqs_to_skip = @_;

    #readin files
    my %full_seq_fastadata = READIN_FASTA($full_seq_fasta);
    my %intein_fastadata = READIN_FASTA($intein_seq_fasta);

    #open output file
    my($file_handle)=($full_seq_fasta=~/.*?\/(.*?)\..*/);
    mkdir($dir_name);
    my $file_to_print = "$dir_name\/$file_handle\_extein.$extension";
    
    open(my $extein_out, "+> $file_to_print");
    
    #for each intein sequence, find the associated full sequence ascession
    foreach my $intein_asc (keys %intein_fastadata){

        #check if current intein was flagged to be skipped
        if(grep(/^$intein_asc$/, @seqs_to_skip)){
            next;
        }

        #find the associated full seq asc
        my $full_asc;
        my $root_key; #not needed for AA
        if($flag == 0){
            #for aa seqs it is:
            $full_asc = $paired_memory{$intein_asc}{"full_aa_asc"};
            if(!defined $full_asc){
                die print "Could not find paired full amino acid accession for $intein_asc!\n\n";
            }
        }
        else{
            #for nucleotides: get key 1
            #print "------------------------------------------------------\n\nNeed to find full nuc asc paired with $intein_asc\n\n";
            foreach my $key1 (keys %paired_memory){
                #print "Checking key $key1 and it's value for intein_nuc_asc $paired_memory{$key1}{\"intein_nuc_asc\"}\n";
                #check if this match was never found
                next if(!defined $paired_memory{$key1}{"intein_nuc_asc"});
                if($paired_memory{$key1}{"intein_nuc_asc"} eq $intein_asc){
                    #print "Found match for $intein_asc!\nAssociated full nuc asc is: $paired_memory{$key1}{\"full_nuc_asc\"}\n\n";
                    $full_asc = $paired_memory{$key1}{"full_nuc_asc"};
                    last;
                }
                else{
                    next;
                }
            }
            if(!defined $full_asc){
                die print "Could not find paired full nucleotide accession for $intein_asc!\n\n";
            }
        }

        #begin to get extein seqs
        my $full_sequence = $full_seq_fastadata{$full_asc};

        #splice the intein out of the full sequence
        my $intein_sequence = $intein_fastadata{$intein_asc};
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
        print $extein_out "$full_asc\n$extein_sequence\n";
    }

    close $extein_out;
}

sub GET_FULL_PROTEIN_SEQ{
    #takes fasta file of intein sequence, and file data hash as input.
    #retrieves the associated complete asc amino acid sequence
    my $fasta = shift;
    my %global_memory_hash;

    #blast file against database
    my($blast_file)=RUN_BLAST("blastp",$prot_database,$fasta,"\"6 qseqid sseqid pident\"");
    
    #READIN query seqs
    my(%fastadata)=READIN_FASTA($fasta);
    
    #parse BLAST output
    my($hashref2,$hashref3)=PARSE_BLAST($blast_file,\%fastadata,0);
    my %blast_matches = %{$hashref2};
    my %local_memory_hash=%{$hashref3};

    #create range file and extract associated sequences
    mkdir("full_prot_seq");
    move($blast_file,"full_prot_seq\/$blast_file");
    my($file_handle)=($fasta=~/.*\/(.*?)\./);
    EXTRACT_MATCHES(\%blast_matches,$prot_database,"full_prot_seq\/$file_handle.faa","0");

    #Now save full protein sequence acessions into the memory hash (dictionary)
    #READIN full sequence fasta file
    my %full_seq_fastadata = READIN_FASTA("full_prot_seq\/$file_handle.faa");
    my @nomatch;
    foreach my $intein_asc (keys %local_memory_hash){
        
        foreach my $full_seq_asc (keys %full_seq_fastadata){

            #get the BLAST full sequence asc 
            my $check = $local_memory_hash{$intein_asc}{"full_seq_blast_asc"};
            
            #compare said asc to full seq asc 
            if($full_seq_asc=~/$check/){
                #save full seq asc and go to next intein asc
                $global_memory_hash{$intein_asc}{"full_aa_asc"}=$full_seq_asc;
                last;
            }
            else{
                next;
            }
        }

        #check if you could not find a match for the intein asc
        if(!defined $global_memory_hash{$intein_asc}{"full_aa_asc"}){
            print "\n\nCould not find match to $intein_asc in full seq aa ascs.\n\n";
            push(@nomatch,$intein_asc);
        }
    }

    return("full_prot_seq\/$file_handle.faa",\%global_memory_hash,@nomatch);
}

sub EXTRACT_MATCHES{
    #takes a hash of matches
    #extracts matches, print them to a range file and extract sequences
    #returns hash of sequences and ascessions
    my (%seqs_for_range)= %{my $hashref = shift};
    my $database = shift;
    my $output_name = shift;
    my $mode = shift;
    my $strand;

    open(my $range, "+> range.txt");
    foreach my $asc (keys %seqs_for_range){
        my @split_data = split(/\t/,$seqs_for_range{$asc});
        
        #if working w/ nucleotides need to include start and end points in range file
        if($mode eq "1"){
            my $plus_or_minus = ($split_data[3]-$split_data[4]);
            if($plus_or_minus >= 0){
                $strand = "minus";
                #print "Line is $_\n";
                print $range "$split_data[1]\ $split_data[4]\-$split_data[3]\ $strand\n";
            }
            else{
                #print "Line is $_\n";
                $strand = "plus";
                print $range "$split_data[1]\ $split_data[3]\-$split_data[4]\ $strand\n";
            }
        }

        #for proteins
        else{
            print $range "$split_data[1]\n";
        }
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
    system("$blastversion -query $query_fasta -db $database -out blast6.txt -evalue 1e-20 -outfmt $out_version");

    return("blast6.txt");
}

sub PARSE_BLAST{
    #INPUTS:
    #blast file outfmt 6
    #fasta sequence data
    #mode variable (0 = protein, 1= nucleotide)

    #OUTPUTD:
    #returns hash of matches (1 per query)
    
    my $blast_to_parse = shift;
    my %sequence_data = %{my $hashref = shift};
    my $mode = shift;

    my $subject_length = 0;
    my %matches;
    my %memory;
    
    #readin BLAST file
    open(my $blast6, "< $blast_to_parse");
    while(<$blast6>){
        chomp;
        my @split = split(/\t/,$_);

        #percent ID cutoff
        if($split[2] <= 99.000){
            #print "Percent ID of match only $split[4] not 100.\n";
            next;
        }

        #check mode
        if($mode == 1){
            $subject_length = abs($split[3]-$split[4])+1; #yes the plus 1 is necessary IDK why.
        }

        #check if match has been previously found to the given query
        if(!defined $matches{$split[0]}){

            #YOU WILL NEED TO EDIT THIS SUBROUTINE TO MATCH YOUR SEQUENCE ACESSION FORMAT
            my($phagename)=GET_ACS_NAME($split[0]);

            if(!defined $phagename){
                die print "Failed to get phage name for $split[0]\n\n";
            }

            #check for matching name
            if($split[1]=~/$phagename/){
                #finds the associated input sequence acession to the current line's subject
                my($paired_intein_match)=FIND_ASSOCIATED_FULL_ASC($split[0],\%sequence_data);

                #check for coverage cutoff if trying to capture nucleotides.
                if($mode == 1){
                    my $query_nucl_length = (length($sequence_data{$paired_intein_match})*3);
                    if($subject_length ne $query_nucl_length){
                        next;
                    }
                    else{
                        $memory{$paired_intein_match}{"query_blast_asc"}=$split[0];
                        $memory{$paired_intein_match}{"subject_blast_asc"}=$split[1];
                        $memory{$paired_intein_match}{"sstart"}=$split[3];
                        $memory{$paired_intein_match}{"send"}=$split[4];
                    }
                }

                #if parsing blastp results, then straight to output
                else{
                    $memory{$paired_intein_match}{"query_blast_asc"}=$split[0];
                    $memory{$paired_intein_match}{"full_seq_blast_asc"}=$split[1];
                    #print "$phagename matches $split[1]!\n";
                }

                #save match
                $matches{$split[0]}=$_;
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
    #takes a string unique to 1 acession line
    #takes a hash with acessions as the keys
    #outputs the acession the string is in

    my $asc_to_find = shift;
    my %ascs = %{my $hashref = shift};
    my $found_it="0";
    
    foreach my $asc_in_file (keys %ascs){
        if($asc_in_file=~/$asc_to_find/){
            $found_it=$asc_in_file;
        }
        else{
            next;
        }
    }

    if($found_it eq "0"){
        $found_it = "NO_MATCH";
    }

    return($found_it);
}

sub GET_ACS_NAME{
    #takes blast query line as input
    #returns the name of the asc or dies
    #You will need to change this subroutine to match the needs of your database

    my $text = shift;
    my $name = "";


    ($name)=($text=~/(.*?)\_.*/);
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

sub APPEND_UID{
    #takes a fasta file
    #appends a unique ID to each acession
    #ALSO CHECKS FOR THE ODD CASES OF A C GETTING PLACED AFTER THE : PRIOR TO THE START AND STOP SITES.
    #NO I DO NOT KNOW WHY BLASTDBCMD DOES THIS FUCKING BS
    my $infile = shift;
    my $unique = 0;
    
    open(IN, "< $infile");
    open(OUT, "+> temp.fasta");
    while(<IN>){
        chomp;
        $_ =~ s/\s+$//;
        if($_=~/\>/){
            #check if the stupid :c is present
            if($_=~/.*?\:c[0-9]+\-[0-9]+/){
                $_ =~ s/\:c/\:/;
            }
            print OUT "$_\_UID\_$unique\n";
            $unique++;
        }
        else{
            print OUT "$_\n";
        }
    }
    close IN;
    close OUT;
    unlink($infile);
    rename("temp.fasta",$infile);
}

sub SETUP{
    #make directories, and get files needed for downstream processes
    my(@intein_fasta)=GET_INPUTS_RECURSIVE($input_directory);

    return(@intein_fasta);

}
