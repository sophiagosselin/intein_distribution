#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;

#need to make phylogenies!
#this is gonna be a process

##First - Need to take extein sequences and run a blast search with them against the phagesdb database
##Each extein should be a separate query

#Second - go through outputs and group together any exteins that match the same results
##This will give us groups of homologous extein sequences to build phylogenies out of

#Third - Import the dictionary file that allows you to match up extein and intein sequences

#Fourth - take intein sequences associated with the exteins in step 1 and cluster them based on identity 

#Fifth - Make a phylogeny out of each group of exteins
##First with the aa files - should be easy to get these sequences (use blastdbcmd plus associated sequences not already in file)
##Second with nuc seqs
###This will require some reusing of code from extract_sequences.pl
###For each protein encoding sequence that does not contain an intein, just use tblastn to find the match
###For those with the intein, use the sequence we have already extracted
###This should help to cut down on the work needed.

#FINALLY - Pair up intein clusts with the associated phylogeny and use the dictionary + R heatmap code to create final figures


### IMPORTANT: MAKE SURE YOUR EXTEIN SEQ FILE HAS THE SAME ACC AS THE ONE IN THE DATABASE. DO NOT MAKE ALTERATIONS.


##GLOBALS
my $extein_aa_file = $ARGV[0]; #these are the extein sequences of proteins known to contain an intein
my $protein_database = $ARGV[1]; #this database needs to have been created with parse_seqids
my $extein_nuc_file = $ARGV[2]; #nuc seqs of the extein seqs of $ARGV[0]
my $dictionary_file = $ARGV[3];

#READIN the dictionary which connects the nuc and aa acc's for inteins and exteins will need for later as global
my(%dictionary)=READIN_DICTIONARY($dictionary_file);


#START
MAIN();

sub MAIN{

    #readin the extein aa file
    my(%extein_aa_data)=READIN_FASTA($extein_aa_file);

    #create array of query names for later parsing
    my(@input_query_accs)=(keys %extein_aa_data);   
    
    #take each extein sequence and output to a seperate file
    my(@extein_aa_files)=SPLIT_FASTA("split_extein_inputs"%extein_aa_data);

    #foreach independent extein fasta file run a blast search and save the file path
    my @extein_aa_blast_files;
    foreach my $extein_aa_query (@extein_aa_files){
        my($blast_file)=BLAST("blastp",$protein_database,$extein_aa_query,"\"6 qseqid sseqid pident qlen slen\"");
        push(@extein_aa_blast_files,$blast_file);
    }

    #readin each BLAST output file and parse the matches 
    #matches that meet the criterion are added to the hash as a subkey
    my %all_extein_aa_blast_results;
    foreach my $blast_file (@extein_aa_blast_files){
        %paired_queries_and_blast_files=PARSE_BLAST($blast_file,\%all_extein_aa_blast_results);
    }

    #go through all BLAST results. Any exteins which have any overlapping matches are grouped together
    my(%grouped_exteins)=FIND_MATCH_OVERLAP(%all_extein_aa_blast_results);

    #take grouped half matrix, and use it to create unified fasta files
    my(@grouped_extein_files)=CREATE_GROUP_FILES("grouped_exteins",\%grouped_exteins);

    #parse each file, and extract matches
    foreach my $group_blast (@grouped_extein_files){
        
        #START - Protein seq processing
        ##remove duplicates and a get a hash of sequences to extract
        my %group_matches;
        (%group_matches)=PARSE_BLAST($group_blast,\%group_matches);

        ##create range file and extract matches from database
        my($main_hashref1,$main_arrayref1)=REMOVE_KEYS(\%group_matches,\@input_query_accs);
        %group_matches = %{$main_hashref1};
        my @seqs_to_add_later = @{$main_arrayref1};

        ##skips matches to query sequences known to contain an intein
        EXTRACT_MATCHES(\%group_matches,"$group_blast.faa");

        #save a copy for alter before adding original extein seqs
        copy("$group_blast.faa","$group_blast.faa.nooriginals");

        ##reunite new extein homologs and original extein sequence
        open(my $reunion, ">> $group_blast.faa");
        foreach my $acc (@seqs_to_add_later){
            print $reunion "$acc\n$extein_aa_data{$acc}\n";
        }
        close $reunion;

        ##MAKE ALIGNMENT
        my($mafft_alignment)=ALIGN_PROT_SEQS("$group_blast.faa");

        ##NOW MAKE A PHYLOGENY HERE
        MAKE_PROT_PHYLOGENY($mafft_alignment);


        #START - Nucleotide seq processing

        #get nucleotide sequences for all proteins in each group
        "$group_blast\.faa"

        #FOR TOMORROW:
        #FINISH THIS SEGMENT OF CODE. $group_blast.faa.nooriginals contains aa seqs that need to be converted into nuc seqs
        #Run a tblastn search with these and extract matches according to the group
        #THEN use the dictionary hash to get the associated nuc seqs for the extein files.
        #when you go to test this code, make sure to remove the "_extein_only" at the end of the accessions
    
    
    }

}

sub MAKE_PROT_PHYLOGENY{
    #takes MAFFT alignment, outputs a phylogenetic tree
    my $inalignment = shift;

    system("iqtree2 -s $inalignment -m MFP -B 1000 -T AUTO");
}


sub ALIGN_PROT_SEQS{
    #takes a multi fasta file
    #makes an alignment using MAFFT

    my $infasta = shift;

    system("mafft --ep 0 --genafpair --maxiterate 1000 $infasta > $infasta\.mafft");

    return("$infasta.mafft");
}

sub EXTRACT_MATCHES{
    #takes a hash of matches
    #extracts matches, print them to a range file and extract sequences
    #returns hash of sequences and ascessions
    my (%seqs_for_range)= %{my $hashref = shift};
    my $output_name = shift;
    my $strand;

    open(my $range, "+> range.txt");
    foreach my $file (keys %seqs_for_range){

        foreach my $match (keys %{$seqs_for_range{$file}}){
            
            print $range "$seqs_for_range{$file}{$match}\n";

        }
    }

    close $range;

    #extract matches
    system("blastdbcmd -db $protein_database  -entry_batch range.txt -outfmt \"%f\" > $output_name");
}

sub REMOVE_KEYS{
    #takes hash and array as inputs.
    #removes keys that match the inputs from hash
    #returns hash

    #imporantly bc this array has subkeys
    my %hash = %{my $hashref1 = shift};
    my @array = @{my $arrayref = shift};
    my @return_array;
    
    foreach my $key_to_remove (@array){

        foreach my $level1 (keys %hash){

            foreach my $level2 (keys %{$hash{$level1}}){

                if($level2 eq $key_to_remove){
                    delete(%{$hash{$level1}{$level2}});
                    push(@return_array,$key_to_remove);
                }
                else{
                    next;
                }
            }
        }
    }

    return(\%hash,\@return_array);
}

sub CREATE_GROUP_FILES{
    #takes half nmatrix as input
    #goes through entry by entry and groups blast outputs together
    #returns directory with subdirectory of new grouped blast files

    my $outdir = shift;
    my %group_data = %{ my $hashref1 = shift};
    my %temporary_groups;
    my $group_counter = 1;

    #go through and find all exteins that should be grouped together for each group
    foreach my $blast_file_1 (keys %group_data){

        foreach my $blast_file_2 (keys %{$group_data{$blast_file_1}}){

            #if 2 exteins belong to the same group
            if($group_data{$blast_file_1}{$blast_file_2} == 1){
                
                #check if one of these files is already in a group
                my $group_found_toggle = 0;
                foreach my $temp_group (keys %temporary_groups){
                    foreach my $index (@{$temporary_groups{$temp_group}}){
                        
                        if($temporary_groups{$temp_group}[$index] eq $blast_file_1){
                            push(@{$temporary_groups{$temp_group}},$blast_file_2);
                            $group_found_toggle = 1;
                            last;
                        }

                        elsif($temporary_groups{$temp_group}[$index] eq $blast_file_2){
                            push(@{$temporary_groups{$temp_group}},$blast_file_1);
                            $group_found_toggle = 1;
                            last;
                        }

                    }

                    if($group_found_toggle == 1){
                        last;
                    }
                }

                #check if one of the files already belonged to a group and the other was assigned
                if($group_found_toggle == 1){
                    next;
                }

                #if not, assign them a new group number
                else{
                    $temporary_groups{$group_counter} => [$blast_file_1,$blast_file_2];
                }

            }

            #they don't belong to the same group
            else{
                next;
            }
            
        }

    }

    #Now that groups have been created, create a unified file for each group
    mkdir($outdir);
    
    #create log file so you can easily track witch blast files belong to which group
    open(my $log, "+> $outdir\/group_log.txt");

    my @group_files;
    foreach my $group (keys %temporary_groups){
        
        #record group #
        mkdir("$oudir\/$group");
        print $log "$group\t";

        #start printing to the new file
        open(my $out, "+> $outdir\/$group\/group_$group\.blast");
        
        foreach my $index (@{$temporary_groups{$group}}){
            #readin file
            my $file = $temporary_groups{$group}[$index];
            #print contents to group file
            open(my $blast, "< $file");
            while(<$blast>){
                print $out $_;
            }
            close $blast;

            #print to log file
            print $log "$blast\t";
        }
        
        close $out;
        push(@group_files,"$outdir\/$group\/group_$group\.blast");
        print $log "\n";
    }

    close $log;

    return(@group_files);
}

sub FIND_MATCH_OVERLAP{
    #takes a hash of queries and matches
    #discover overlap between queries, and return the set of blast files to be merged into clusters/groups
    my %matches = shift;

    #data structure for each comparison should be as follows:
    #hash{obj1}{obj2} = BOOLEAN (true for same group, false for different)
    #can avoid extra work by checking if hash{obj1}{obj2} exists when starting with obj2 as the starting point
    my %grouping_hash;

    #the first level key is the blast file
    foreach my $blast_file_1 (keys %matches){

        #compare each one to all other blast files
        foreach my $blast_file_2 (keys %matches){

            #skip self-self comparisons
            if($blast_file_1 eq $blast_file_2){
                next;
            }

            #skip redundant work
            if(exists $grouping_hash{$blast_file_2}{$blast_file_1}){
                next;
            }

            #toggle if a match is found
            $toggle = 0;

            #now check if there are any overlapping matches between the 2 files
            foreach my $blast_match_1 (keys %{$matches{$blast_file_1}}){
                foreach my $blast_match_2 (keys %{$matches{$blat_file_2}}){
                    
                    #there is overlap! Record match
                    if($blast_match_1 eq $blast_match_2){
                        $grouping_hash{$blast_file_1}{$blast_file_2} = 1;                      
                        $toggle = 1;
                        last;
                    }

                    else{
                        next;
                    }
                }
            }

            #if there are no matches, record that!
            if($toggle == 0){
                $grouping_hash{$blast_file_1}{$blast_file_2} = 0;
            }
        }
    }

    return(%grouping_hash);
}

sub PARSE_BLAST{
    #takes a blast file outfmt 6
    #returns hash of matches (1 per subject ID)

    my $blast_to_parse = shift;
    my %matches = %{my $hashref = shift};
    my $subject_length = 0;
    
    #readin BLAST file
    open(my $blast6, "< $blast_to_parse");
    while(<$blast6>){
        chomp;
        my @split = split(/\t/,$_);

        #check if this is the first match by this subject hit
        if(exists $matches{"$blast_to_parse"}{$split[1]}){
            next;
        }
        
        #check if it meets the coverage criterion
        my $percent_cutoff = ($split[3]*0.5);
        if($split[4]>=$percent_cutoff){
            
            #if it does, record match!
            $matches{"$blast_to_parse"}{$split[1]}=$split[0];

        }

        #otherwise next entry!
        else{
            next;
        }

    }

    close $blast6;

    return(%matches,);
}

sub BLAST{
    #take a blastversion (command), database, query, and out fmt and runs BLAST search 
    my ($blastversion) = shift;
    my ($database) = shift;
    my ($query_fasta) = shift;
    my ($out_version) = shift;
    system("$blastversion -query $query_fasta -db $database -out $query_fasta.blast -evalue 1e-30 -outfmt $out_version");

    return("$query_fasta.blast");

}

sub SPLIT_FASTA{
    #takes the name of an output directory
    #takes a hash of accessions and sequences
    #prints each one to a file
    my $outdir = shift;
    my %seq_data = %{my $hashref1 = shift};
    my @files;

    mkdir($outdir);

    #print each seq and accession to an independent file
    foreach my $acc (keys %seq_data){
        open(my $out, "+> $outdir\/$acc\.fasta");
        print my $out "$acc\n$seq_data{$acc}\n";
        close $out;
        #save to array
        push(@files,"$outdir\/$acc\.fasta");
    }

    return(@files);
}


sub READIN_FASTA{
    #takes fasta as input
    #returns hash using accessions for keys and sequences for values 
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