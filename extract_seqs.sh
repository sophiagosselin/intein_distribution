#!/usr/bin/env bash
#SBATCH --job-name=phagesdb_get_intein_ass0ciated_seqs
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
export BLASTDB=$BLASTDB:/home/FCAM/sgosselin/phamdb/


perl extract_sequences.pl intein_prot_seq_clusters/ allphams.faa phagesdb.fna nt