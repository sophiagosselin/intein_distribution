#!/usr/bin/env bash
#SBATCH --job-name=intein_geo_dist
#SBATCH --nodes=1
#SBATCH --qos=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32gb
#SBATCH --mail-type=END
#SBATCH --mail-user=sophia.gosselin@uconn.edu
#SBATCH -o geo_%j.out
#SBATCH -e geo_%j.err

echo hostname

#dependencies
module load blast/2.11.0
module load muscle
module load perl/5.36.0
module load usearch

#add libraries to path
export PERL5LIB=/home/FCAM/sgosselin/perl5/lib/perl5
export export PATH=~/home/FCAM/sgosselin/phamdb:$PATH


perl geo_dist_pipe.pl all_matches.fasta .975 32 /home/FCAM/sgosselin/phamdb/allphams.faa
