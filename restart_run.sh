#!/bin/bash
###Define paths
#run_path= pwd

###Go to run directory
#cd run_path


comm1=`sbatch runSCHISM.sh /scratch/g/g260204/test/test/sgtest014_2010`
IFS=' ' read -ra ADDR <<< "$comm1"  #Internal Field Seperator#
id=${ADDR[3]}

#id=19707501

#for n in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24; do
for n in {1..52} ; do
  #comm2=`sbatch merge_hotstart.sh`
  comm2=`sbatch --depend=afterany:$id merge_hotstart.sh $n`
  IFS=' ' read -ra ADDR <<< "$comm2"  #Internal Field Seperator
  id=${ADDR[3]}

#  comm3=`sbatch --depend=afterany:$id copy_mirror.sh $n`
#  IFS=' ' read -ra ADDR <<< "$comm3"  #Internal Field Seperator
#  id=${ADDR[3]}

  comm4=`sbatch --depend=afterany:$id runSCHISM_hot.sh /scratch/g/g260204/test/test/sgtest014_2010 $n`
  IFS=' ' read -ra ADDR <<< "$comm4"  #Internal Field Seperator
  id=${ADDR[3]}
done

