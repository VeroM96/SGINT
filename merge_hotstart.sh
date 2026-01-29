#!/bin/bash
#SBATCH --job-name=mhot_sgtest014_2010
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:40:00
#SBATCH --mail-type=NONE
#SBATCH --account=gg1302
#SBATCH --output=merge_hotstart.o
#SBATCH --error=merge_hotstart.e

###Load  module files
module load netcdf_c/4.3.2-gcc48
module load intel/18.0.4
module load intelmpi/5.1.3.223
module load python3
export I_MPI_PMI_LIBRARY=/use/lib64/libmpi.so
export LD_LIBRARY_PATH="/sw/spack-levante/mambaforge-22.9.0-2-Linux-x86_64-kptncg/lib:$LD_LIBRARY_PATH"

#read arguments
ver=$1 #SGINT step; starting from 1 after 
#check for arguments, else stop script
if [ -z "$ver" ]
then
  echo "No SGINT step argument supplied, exiting"
  exit 1
fi

#check if seagrass files are in folder
if [ ! -f sav_N.gr3 ] || [ ! -f sav_RB.gr3 ] || [ ! -f sav_LB.gr3 ] || [ ! -f sav_h.gr3 ] 
then
  echo "Seagrass files for step $ver not found, exiting"
  exit 1

#define important directories
this_dir=$PWD
merge_bin_path=/work/gg0877/g260204/schism_ecosmo_sg/schism_v5.11.0/build/bin/combine_hotstart7
mergepath=/scratch/g/g260204/test/test/sgtest014_2010
cd $mergepath

#find biggest timestep
biggest_t=0
for f in hotstart_*_*.nc; do
  file=$f
  IFS='_' read -ra ADDR <<< "$file"  #Internal Field Seperator
  dtnc=${ADDR[2]}
  IFS='.' read -ra ADDR <<< "$dtnc"  #Internal Field Seperator
  dt=${ADDR[0]}
  if [ $biggest_t -lt $dt ]
    then biggest_t=$dt
  fi
done

echo "Biggest timestep found: $biggest_t"

#merge hotstart with biggest timestep and copy it here
$merge_bin_path -i $biggest_t
rm $this_dir\/hotstart.nc
ln -s $this_dir\/outputs\/hotstart_it\=$biggest_t\.nc $this_dir\/hotstart.nc
cd $this_dir
echo "Hotstart merged and linked to $this_dir/hotstart.nc"

#copy mirror file
new="outputs/mirror_${ver}.out"
cp outputs/mirror.out $new
echo "Mirror file copied to $new"

#calculate means for sgint 
sgintmeanpath=/home/g/g260204/schism_v5.11/sgint_scripts_clean/make_week_mean_parallel.py
python3 $sgintmeanpath $this_dir $ver
echo "Week mean calculated for $ver"

#calculate sg int
sgintpath=/home/g/g260204/schism_v5.11/sgint_scripts_clean/run_sgint_parallel.py
python3 $sgintpath $this_dir $ver
echo "SG int calculated for $ver"

#calculate POC source
makesourcepath=/home/g/g260204/schism_v5.11/sgint_scripts_clean/create_source/create_source_nc_sep.py
python3 $makesourcepath $this_dir $ver
rm $this_dir\/source.nc
ln -s $this_dir\/outputs/week_mean/source_$ver.nc $this_dir\/source.nc
echo "source file created for $ver"

#calculate POC sediment changes
cp $this_dir\/outputs\/hotstart_it\=$biggest_t\.nc $this_dir\/outputs\/hotstart_it\=$biggest_t\_orig.nc
updatehotpath=/home/g/g260204/schism_v5.11/sgint_scripts_clean/create_sed_source/sed_source_sep.py
python3 $updatehotpath $this_dir $ver
echo "updated sediment fraction in hotstart"
