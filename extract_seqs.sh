#!/usr/bin/env bash
#SBATCH --job-name=peanuts_get_seq_nolim
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o extract_%j.out
#SBATCH -e extract_%j.err

echo hostname

#dependencies
module load blast/2.11.0
module load perl/5.36.0

#add libraries to path
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/phagesdb/
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/allphams/


perl extract_sequences_org.pl inputs/ /home/FCAM/sgosselin/allphams/allphams.faa /home/FCAM/sgosselin/phagesdb/phagesdb.fna
