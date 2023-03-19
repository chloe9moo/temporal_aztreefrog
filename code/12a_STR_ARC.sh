#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=4
#SBATCH -t 12:00:00
#SBATCH -J str_run
#SBATCH -p normal_q
#SBATCH --account=usgs_rcs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=chloe9mo@vt.edu

#TITLE: STRUCTURE
#AUTHOR: C. E. Moore
#UPDATED: 14 NOV 2022

#AIM: Iterate through values of K, with multiple replicates
#set up run (not on cluster only)
#cd ~/Documents/Projects/console
#rm structure_loci.txt
cp ../microsat_data/structure_loci_45678.txt ../../console/
cd ~/Documents/Projects/console
sed -i '1d' structure_loci_45678.txt #remove locus names
#variables
NUMK=6 #for full run, number of pops + 1
REPS=10 #for full run, 10 reps of each K
STR_FI='structure_loci_45678.txt'
NUMIND=346

#make function for loop
iter_struct () {
	local k=$1
	for i in $(seq 1 $REPS)
	do
		./structure > runK${k}_r$i.txt -m mainparams -K $k -N $NUMIND -i $STR_FI -o Year_Pop_K${k}_r$i
        	mv Year_Pop_K${k}_r${i}_f Output_Files/
		mv runK${k}_r$i.txt Run_Files/
	done
	wait
        echo "STRUCTURE with K=$k complete"
}

#run program
for k in $(seq 1 $NUMK)
do
	iter_struct "$k" &
done
wait

echo "STRUCTURE runs complete..."

#pull in results
#rm iter_structure_results.csv
printf "K\tRep\tMeanLnP_K\n" >> iter_structure_results_45678.csv #make csv

for k in $(seq 1 $NUMK)
do
	for rep in $(seq 1 $REPS)
	do
		#get line with LnP(K) value
		LINE="$(grep -hnr "Estimated Ln Prob of Data" Output_Files/Year_Pop_K${k}_r${rep}_f)"
		#get value
		LnK=${LINE#*= }
		#add to csv
		printf "$k\t$rep\t$LnK\n" >> iter_structure_results_45678.csv
	done
	wait
	
	echo "K=$k added to csv..."
done
wait

echo "Results csv complete..."




