#!/bin/bash
#
#SBATCH --job-name=Folding
#SBATCH --partition=compute
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --output=TranscriptFold.j%.out
#SBATCH --error=TranscriptFold.j%.err
#SBATCH --mail-type=all
#SBATCH --mail-user=lelal@clemson.edu
#
#Add required module for RNA folding and path to data tables
module add rnastructure

export DATAPATH=/opt/ohpc/pub/Software/RNAstructure/data_tables
export OMP_NUM_THREADS=10

#run rnastructure for every seven lines
while read line
  do
  #echo $line
  arr=($line)
  header=$(printf '>%s' "$line")
  name=${arr[0]}
  sequence=${arr[1]}
  echo -e ">${name}" > temp${name}.fasta
  echo $sequence >> temp${name}.fasta
  echo ${name}.start
  Fold-smp temp${name}.fasta ${name}.ct -mfe -md 200
  partition-smp temp${name}.fasta ${name}.pfs -md 200
  #ProbabilityPlot temp${name}.pfs temp${name}.dp -t
  MaxExpect ${name}.pfs ${name}_prob.ct
  efn2 ${name}.ct ${name}_mfe_energylist.txt
  efn2 ${name}_prob.ct ${name}_prob_energylist.txt
  echo ${name}.done
done < 9-10-21_TranscriptsFold.tab

module unload rnastructure
