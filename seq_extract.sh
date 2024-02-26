#!/usr/bin/env bash
#SBATCH --job-name=extract_seqs
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16gb
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o seqs_%j.out
#SBATCH -e seqs_%j.err

echo hostname

#dependencies
module load blast/2.11.0
module load perl

#add libraries to path
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/allphams/
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/phagesdb/

perl extract_sequences.pl inputs/ allphams.faa phagesdb.fna
