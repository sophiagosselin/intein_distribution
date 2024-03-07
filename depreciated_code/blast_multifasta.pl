#!/usr/bin/perl -w
use strict;
use warnings;

#takes input multifasta breaks it up and runs individual blast searches with each
#captures homolgous sequences for the purpose of alignemt investigation

my $input_fasta_file = $ARGV[0];

my $database = $ARGV[1];

MAIN();

sub MAIN{
    my(%file_and_phams) = BREAKUP_MULTIFASTA($input_fasta_file);

    foreach my $file (keys %file_and_phams){
        mkdir("blast_files");
        my ($file_handle)=($file=~/.*?\/(.*?)\.fasta/);
        my $output_blast = RUN_BLAST($file,$file_handle);
        mkdir("extracted_matches");
        EXTRACT_MATCHES($output_blast,$file_handle);
    }
}

sub EXTRACT_MATCHES{
    my $blast_results = shift;
    my $output_path = shift;
    my %matches;

    open(my $range, "+> range.txt");
    open(my $blast, "< $blast_results");
    while(<$blast>){
        chomp;
        my @split_lines = split(/\t/,$_);
        if(defined $matches{$split_lines[1]}){
            next;
        }
        else{
            print $range "$split_lines[1]\n";
            $matches{$split_lines[1]}=1;
        }
    }
    close $blast;
    close $range;

    system("blastdbcmd -db $database -entry_batch range.txt -outfmt \"%f\" > extracted_matches\/$output_path\.extracted\.fasta");
}

sub RUN_BLAST{
    my $query_fasta = shift;
    my $handle = shift;
    my $out_name = "blast_files\/$handle\_blast6.txt";
    system("blastp -query $query_fasta -db $database -out $out_name -evalue 1e-10 -outfmt \"6 qseqid sseqid pident\"");
    return($out_name);
}

sub BREAKUP_MULTIFASTA{
    my $input_file = shift;
    my %paired_files_and_phams;
    mkdir("single_fastas");


    my %multifasta_data = READIN_FASTA($input_file);
    foreach my $sequence_asc (keys %multifasta_data){
        my($seq_name,$pham)=($sequence_asc=~/\>(.*?)\ .*?\ (.*?)\ .*/);        
        open(my $single_fasta, "+> single_fastas\/$seq_name\.fasta");
        print $single_fasta "$sequence_asc\n$multifasta_data{$sequence_asc}\n";
        close $single_fasta;
        $paired_files_and_phams{"single_fastas\/$seq_name\.fasta"}=$pham;
    }

    return(%paired_files_and_phams);
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