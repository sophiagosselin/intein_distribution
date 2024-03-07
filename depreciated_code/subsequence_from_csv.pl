#!/usr/bin/perl -w
use strict;
use warnings;

#takes file with list of phage names, and start/stop points.
#finds genome in database, and extracts seq of interest

my $database = $ARGV[0];
my $input_csv = $ARGV[1];


MAIN();

sub MAIN{
    #open input file and get seqs of interest.
    #make sure there isn't a header, and the csv has the following format:
    #phage name \t Seq Start \t Seq  End
    open(my $csv, "< $input_csv");
    open(my $out, "+> sequences_from_db.fasta");
    while(<$csv>){
        chomp;
        my $seq_of_interest = "";
        my @data = split(/\t/,$_);
        my $phage_name = $data[0];
        my $seq_start = $data[1];
        my $seq_end = $data[2];
        my ($ascession,$genome) = FIND_PHAGE($phage_name);
        if($ascession eq "NA"){
            next;
        }
        my $length = $seq_start-$seq_end; 
        if($length <= 0){
            #REVERSE IT
            $length = abs($length);
            my $seq_wrong_strand = substr($genome,$seq_start,$length) or die print "FUCK YOU BELFORT\nStart: $seq_start\nStop: $seq_end\n$genome\n\n";
            my @nucleotides = split(//,$seq_wrong_strand);
            foreach my $nuc (@nucleotides){
                $nuc = uc($nuc);
                if($nuc eq "A"){
                    $seq_of_interest="A".$seq_of_interest;
                }
                elsif($nuc eq "T"){
                    $seq_of_interest="T".$seq_of_interest;
                }
                elsif($nuc eq "G"){
                    $seq_of_interest="G".$seq_of_interest;
                }
                elsif($nuc eq "C"){
                    $seq_of_interest="C".$seq_of_interest;
                }
                else{
                    die print "Unkown nucleotide! - $nuc\nFrom asc $ascession\nLenght $length\nStart: $seq_start\tEnd:$seq_end\nSequence:\n$seq_wrong_strand\nGenome is\n$genome\n\n";
                }
            }
        }
        else{
            $seq_of_interest = substr($genome,$seq_start,$length);
        }
        print $out "$ascession\t$seq_start\t$seq_end\n$seq_of_interest\n";
    }

    close $csv;
    close $out;
}



sub FIND_PHAGE{
    #find genome in database
    my $name_to_find = shift;
    my $toggle = 0;
    my ($asc,$gen) = "";

    open(my $db, "< $database");
    while(<$db>){
        $_=~s/\s//g;
        chomp;
        if($_=~/\>/){

            my $test_asc = uc($_);
            my $uc_name_to_find = uc($name_to_find);
            if($test_asc=~/\>$uc_name_to_find?\_.*/){
                $toggle = 1;
                $asc = $_;
            }

            elsif($toggle == 1){
                last;
            }

            else{
                next;
            }
        }
        else{

            if($toggle == 0){
                next;
            }

            else{
                $gen.= $_;
            }
        }
    }
    close $db;

    if($gen eq ""){
        print "Could not extract genome for $name_to_find!\n";
        $asc = "NA";
    }
    return($asc,$gen);
}