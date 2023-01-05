#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=40G
#SBATCH -t 48:00:00
#SBATCH -J pred_run
#SBATCH -p normal_q
#SBATCH --account=usgs_rcs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=chloe9mo@vt.edu

cd /fastscratch/chloe9mo/temporal_aztreefrog

module load containers/singularity
echo "Modules loaded:"
module list
echo " "
echo "============================="
echo "Running from:"
pwd
echo " "
echo "============================="


echo "Running AMOVA..."
echo "============================="
singularity exec --bind=/fastscratch/chloe9mo/temporal_aztreefrog:/data /projects/arcsingularity/ood-rstudio141717-geospatial_4.1.1.sif Rscript /data/11a_AMOVA_ARC.R
echo "============================="
echo " "




