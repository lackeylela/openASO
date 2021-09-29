#!/bin/bash
#
#SBATCH --job-name=Folding
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --output=out.%j%.out
#SBATCH --error=error.%j%.err
#SBATCH --mail-type=all
#SBATCH --mail-user=
#
#Add required module for RNA folding and path to data tables
module add rnastructure

export DATAPATH=/opt/ohpc/pub/Software/RNAstructure/data_tables

for file in *ct; do
  echo $file
  efn2 $file ${file}_mod_energylist.txt
done

module unload rnastructure
