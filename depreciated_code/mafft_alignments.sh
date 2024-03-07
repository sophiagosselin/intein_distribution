#!/usr/bin/env bash
#SBATCH --job-name=mafft_alignments
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4gb
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o mafft_%j.out
#SBATCH -e mafft_%j.err

echo hostname
module load mafft

sequence_files=(extracted_matches/*)

#create directory for ICE BLAST and prep variables
mkdir "MAFFT_alignments"

for file in ${sequence_files[@]}
do
    [[ $file =~ .*\/(.*)\.fasta ]]
    file_handle=${BASH_REMATCH[1]}
    mafft --ep 0 --genafpair --maxiterate 1000 $file > MAFFT_alignments/$file_handle.mafft
done