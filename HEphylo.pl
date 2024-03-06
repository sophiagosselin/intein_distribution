#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;

# HEphylo.pl - (H)omologous (E)xtein phylogeny maker
#02_28_24 By Sophia Gosselin
# takes a nucleotide and amino acid database and a dataset of intein containing exteins
# the extein files must already have the inteins removed
# this script takes those files, finds homologous sequences within the database to act as outgroups, and creates phylogenies of extein seqs.

# IF YOU ARE RUNNING THIS, HEAD THE FOLLOWING WARNINGS:
## IMPORTANT: MAKE SURE YOUR EXTEIN SEQ FILE HAS THE SAME ACC AS THE ONE IN THE DATABASE. DO NOT MAKE ALTERATIONS.
## YOU WILL NEED TO EDIT SOME LINES TO MATCH THE FORMAT OF YOUR ACCESSIONS. SPECIFICALLY: 276
## READ THE DESCRIPTIONS FOR THE INPUTS BELOW. SOME HAVE SPECIFIC REQS.

##GLOBALS
my $extein_aa_file = $ARGV[0]; #these are the extein sequences of proteins known to contain an intein
my $protein_database = $ARGV[1]; #this database needs to have been created with parse_seqids
my $extein_nuc_file = $ARGV[2]; #nuc seqs of the extein seqs of $ARGV[0]
my $nucleotide_database = $ARGV[3]; #this database needs to have been created with parse_seqids
my $paired_dictionary_file = $ARGV[4]; #this file must be in "extien_aa_acc  exteintein_nuc_acc\n" format such that we can link the two input files
my %query_dictionary; 
#this hash will save you a lot of hastle.
    #formatted as follows:
    #Level 1 key: the extein aa original asc
    #Level 2 keys:
        #1. associated split fasta file
        #2. associated blastp file for said split file
        #3. its assigned group number

#START
MAIN();

sub MAIN{

    #readin the extein aa and nuc files
    my(%extein_aa_data)=READIN_FASTA($extein_aa_file);
    my(%extein_nuc_data)=READIN_FASTA($extein_nuc_file);

    #create array of query names for later parsing
    my(@input_aa_accs)=(keys %extein_aa_data);    
    
    #take each extein sequence and output to a seperate file - file names are inside of dictionary file.
    SPLIT_FASTA("split_extein_inputs",\%extein_aa_data);

    #foreach independent extein fasta file do the following
    my %all_extein_aa_blast_results;
    foreach my $extein_aa_query (keys %query_dictionary){

        #run a blast search and save the file path
        my($blast_file)=BLAST("blastp",$protein_database,$query_dictionary{"$extein_aa_query"}{"split_file"},"\"6 qseqid sseqid pident qlen slen\"");
        $query_dictionary{$extein_aa_query}{"blastp"}=$blast_file;

        #readin each BLAST output file and parse the matches 
        #matches that meet the criterion are added to the hash as a subkey
        (%all_extein_aa_blast_results)=PARSE_BLAST($blast_file,\%all_extein_aa_blast_results,0);
    }

    #go through all BLAST results. Any exteins which have any overlapping matches are grouped together
    my(%grouped_exteins)=FIND_MATCH_OVERLAP(\%all_extein_aa_blast_results);

    #take grouped half matrix, and use it to create unified fasta files
    my(@grouped_extein_blastfile)=CREATE_GROUP_FILES("grouped_exteins",\%grouped_exteins);

    #parse each file, and extract matches
    foreach my $group_blast (@grouped_extein_blastfile){
        
        #################################
        #START - Protein seq processing

        ##remove duplicates and a get a hash of sequences to extract
        my %dummyhash;
        my(%group_matches)=PARSE_BLAST($group_blast,\%dummyhash,0);

        ##skips matches to query sequences known to contain an intein
        my($main_hashref1)=REMOVE_KEYS(\%group_matches,\@input_aa_accs);
        %group_matches = %{$main_hashref1};

        ##create range file and extract matches from database
        EXTRACT_MATCHES(\%group_matches,"$group_blast.faa",$protein_database,0);

        #### MIGHT BE A GOOD IDEA TO CLUSTER THESE SEQUENCES AT THIS JUNCTURE.
        my($clustered)=CLUSTER_SEQS("$group_blast.faa");

        #save a copy for later before adding original extein seqs
        copy("$clustered","$group_blast.faa.nooriginals");

        ##reunite new extein homologs and original extein sequence
        open(my $reunion, ">> $clustered");
        my($group_number)=($group_blast=~/.*?\/group_(.*?)\.blast/);
        foreach my $extein_aa_acc (keys %query_dictionary){
            my $group_to_check = $query_dictionary{$extein_aa_acc}{"group_number"};
            #print "\nDebug test:\nGroup number is $group_number from $group_blast!\nValue to compare to is $extein_group_assignments{$extein_in_group}\nKey for that value was $extein_in_group\n\n";
            if($group_number eq $group_to_check){
                
                print $reunion "$extein_aa_acc\n$extein_aa_data{$extein_aa_acc}\n";
            }
            else{
                next;
            }
            
        }
        close $reunion;

        ##MAKE ALIGNMENT
        my($mafft_alignment)=ALIGN_PROT_SEQS("$clustered");

        ##NOW MAKE A PHYLOGENY HERE
        MAKE_PROT_PHYLOGENY($mafft_alignment);

        ##################################
        #START - Nucleotide seq processing

        #get nucleotide sequences for all proteins in each group
        my($tblastn_file)=BLAST("tblastn",$nucleotide_database,"$group_blast.faa.nooriginals","\"6 qseqid sseqid pident qlen sstart send\"");

        ##get best hit for each seq.
        my %group_tblastn_matches;
        (%group_tblastn_matches)=PARSE_BLAST($tblastn_file,\%group_tblastn_matches,1);

        ##create range file and extract matches from database
        EXTRACT_MATCHES(\%group_tblastn_matches,"$group_blast.fna",$nucleotide_database,1);

        ##reunite extein sequences using the dictionary file
        my(%paired_aa_nuc)=READIN_DICTIONARY($paired_dictionary_file);

        open(my $nuc_reunion, ">> $group_blast.fna");
        foreach my $extein_aa_acc (keys %query_dictionary){
            my $group_to_check = $query_dictionary{$extein_aa_acc}{"group_number"};
            if($group_number eq $group_to_check){
                my $nuc_acc = $paired_aa_nuc{"$extein_aa_acc"};
                my $nuc_seq = $extein_nuc_data{"$nuc_acc"};
                print $nuc_reunion "$nuc_acc\n$nuc_seq\n";
            }
            else{
                next;
            }
            
        }
        close $nuc_reunion;

        ##MAKE ALIGNMENT
        my($nuc_mafft_alignment)=ALIGN_PROT_SEQS("$group_blast.fna");

        ##NOW MAKE A PHYLOGENY HERE
        MAKE_PROT_PHYLOGENY($nuc_mafft_alignment);
    }

}

sub INPUT_DICTIONARY{
    #takes 3 inputs
    #current dictionary hash
    #new hash of lvl 1 keys and new values
    #the name of the lvl 2 key to place the new values under
    my %current_dic = %{my $hashref1 = shift};
    my %new_pairs = %{my $hashref2 = shift};
    my $new_key = shift;

    foreach my $level1_key (keys %new_pairs){
        $current_dic{$level1_key}{"$new_key"}=$new_pairs{$level1_key};
        #print "\n\nDebugging: New Dictionary Entry.\nLevel 1 key: $level1_key\nLevel 2 key: $new_key\nLevel 2 value: $new_pairs{$level1_key}\n";
    }
    
    return(%current_dic);
}

sub READIN_DICTIONARY{
    #this sub assumes the aa sequence acc is the first entry and the nucleotide is the second.
    my $infile = shift;
    my %dictionary;

    open(my $dic, "< $infile");
    while(<$dic>){
        if($_=~/\>/){
            my @entries = split(/\t/,$_);
            chomp($entries[1]);
            $dictionary{"$entries[0]"}="$entries[1]";
        }
        else{
            #this line is likely a header or something else
            next;
        }
    }
    close $dic;

    return(%dictionary);
}

sub CLUSTER_SEQS{
    #takes fasta file and runs through usearch at 97.5% ID to cluster out very closely related sequences
    my $infasta = shift;

    system("usearch --cluster_fast $infasta --sort length --id 0.975 --centroids $infasta.clustered");

    return("$infasta.clustered");
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
    my $extractdb = shift;
    my $mode = shift;
    my $strand;

    open(my $range, "+> $output_name.range");
    foreach my $file (keys %seqs_for_range){

        foreach my $match (keys %{$seqs_for_range{$file}}){
            
            #aa mode
            if($mode == 0){
                print $range "$match\n";
            }

            #nuc mode
            elsif($mode ==1){
                my $toprint = $seqs_for_range{"$file"}{"$match"};
                print $range "$toprint";
            }

            else{
                die print "You didn't supply a mode??? - Please contact developer.\n";
            }

        }
    }

    close $range;

    #extract matches
    system("blastdbcmd -db $extractdb -entry_batch $output_name.range -outfmt \"%f\" > $output_name");
}

sub REMOVE_KEYS{
    #takes hash and array as inputs.
    #removes keys that match the inputs from hash
    #returns hash

    #imporantly bc this array has subkeys
    my %hash = %{my $hashref1 = shift};
    my @array = @{my $arrayref = shift};
    
    foreach my $key_to_remove (@array){

        foreach my $level1 (keys %hash){

            foreach my $level2 (keys %{$hash{$level1}}){
                
                #this regex needs to take the original acc saved as key_to_remove and format it in such that it matches the call BLAST uses for it
                #i.e. Original acc: ">Yourseqname proteinannotation somenumbermaybe" BLAST acc: "Yourseqname"
                #Hence you need to make a regex that takes the first and turns it into the second, such that we can properly search for the query. 
                my($comparator)=($key_to_remove=~/\>(.*?)\ .*/);
                #print "\nComparing $level2 to $key_to_remove\nSpecifically the shortened $comparator\n";

                if($level2 eq $comparator){
                    print "They are equal!!\n";
                    delete($hash{$level1}{$level2});
                }
                else{
                    #print "They are not equal :(\n";
                    next;
                }
            }
        }
    }

    return(\%hash);
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
        my $didwefindagroup = 0;

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
                    $didwefindagroup = 1;
                    next;
                }

                #if not, assign them a new group number
                else{
                    my @temp = ($blast_file_1,$blast_file_2);
                    $temporary_groups{$group_counter} = [ @temp ];
                    $group_counter++;
                }

            }

            #they don't belong to the same group
            else{
                next;
            }
            
        }

        #this checks if blast_file_1 is a singleton or not
        if($didwefindagroup == 1){
            next;
        }

        else{
            #it is a singleton! assign a new group!
            $temporary_groups{$group_counter} = [$blast_file_1];
            $group_counter++;
        }

    }

    #Now that groups have been created, create a unified file for each group
    mkdir($outdir);
    
    #create log file so you can easily track witch blast files belong to which group
    open(my $log, "+> $outdir\/group_log.txt");

    my @group_files;
    foreach my $group (keys %temporary_groups){
        
        #record group #
        mkdir("$outdir\/$group");
        print $log "$group\t";

        #start printing to the new file
        open(my $out, "+> $outdir\/$group\/group_$group\.blast");
        
        foreach my $index (@{$temporary_groups{$group}}){
            #readin file
            #print contents to group file
            open(my $blast, "< $index");
            while(<$blast>){
                print $out $_;
            }
            close $blast;

            #print to log file
            print $log "$blast\t";
            
            #save this group number
            #need the regex to recover original seq acc
            foreach my $topquery (keys %query_dictionary){
                if($index eq $query_dictionary{$topquery}{"blastp"}){
                    $query_dictionary{$topquery}{"group_number"}=$group;
                }

            }
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
    my %matches = %{my $hashref = shift};

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
            my $toggle = 0;

            #now check if there are any overlapping matches between the 2 files
            foreach my $blast_match_1 (keys %{$matches{$blast_file_1}}){
                foreach my $blast_match_2 (keys %{$matches{$blast_file_2}}){
                    
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
    my $mode = shift; #this just adds a seperate check for nucleotide seqs to ensure that 1. Only the best sequence if being captured, 2. That only 1 match per subject is captured, and 3. that the coverage is not 150-50% but nearly 100%
    my $subject_length = 0;
    my %nuc_perc_id;
    
    #readin BLAST file
    open(my $blast6, "< $blast_to_parse");
    while(<$blast6>){
        chomp;
        my @split = split(/\t/,$_);
        
        #amino acid results mode
        if($mode == 0){
            #check if this is the first match by this subject hit
            if(exists $matches{"$blast_to_parse"}{"$split[1]"}){
                next;
            }

            #check if it meets the coverage criterion
            my $coverage_cutoff_low = ($split[3]*0.5);
            my $coverage_cutoff_high = ($split[3]*1.5);

            #if it does, record match!
            if($split[4]>=$coverage_cutoff_low && $split[4]<=$coverage_cutoff_high){
                $matches{"$blast_to_parse"}{"$split[1]"}=$split[0];
            }
            
            #otherwise next entry!
            else{
                next;
            }

        }

        #nucleotide results mode
        else{
            #get coverage value
            my $strand = "minus";
            my ($start,$end);
            my $length = $split[4]-$split[5];
            if($length >= 0){
                #$strand = "minus";
                $start = $split[5];
                $end = $split[4];
                
            }
            else{
                $strand = "plus";
                $length = abs($length);
                $start = $split[4];
                $end = $split[5];
            }

            my $cv = ($length/($split[3]*3));
            #check for previous matches with better results
            if(exists $nuc_perc_id{"$split[0]"}){
                #check if this sequence has already been matched with a higher percent ID hit
                if($nuc_perc_id{"$split[0]"}{"pid"} > $split[2]){
                    next;
                }
                #check if this sequences has already been matched with a higher percent cv hit
                if($nuc_perc_id{"$split[0]"}{"cv"} > $cv){
                    next;
                }
            }

            #check if it meets the percent id criterion
            if($split[2] <= 99.000){
                next;
            }

            #check coverage cutoff and record match if it makes the cut.
            if($cv >= 0.99000 && $cv <= 1.01000){
                $matches{"$blast_to_parse"}{"$split[0]"}="$split[1]\ $start\-$end\ $strand\n";
                $nuc_perc_id{"$split[0]"}{"pid"} = $split[2];
                $nuc_perc_id{"$split[0]"}{"cv"} = $cv;
            }

            else{
                next;
            }

        }

    }

    close $blast6;

    return(%matches);
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
    my %paired_files;

    mkdir($outdir);

    #print each seq and accession to an independent file
    foreach my $acc (keys %seq_data){
        my($file_name)=($acc=~/\>(.*)/);
        $file_name=~s/\ /\_/g;

        open(my $split, "+> $outdir\/$file_name\.fasta");
        print $split "$acc\n$seq_data{$acc}\n";
        close $split;
        #save to array
        $paired_files{$acc}="$outdir\/$file_name\.fasta";
    }

    #add to dictionary!
    %query_dictionary = INPUT_DICTIONARY(\%query_dictionary,\%paired_files,"split_file");
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
        if($_=~/\>/){
            $accession=$_;
            $sequences{"$accession"}="";
        }
        else{
            $sequences{"$accession"}.=$_;
        }
    }
    close IN;
    return(%sequences);
}

