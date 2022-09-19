#!/bin/bash
#TITLE: STRUCTURE
#AUTHOR: C. E. Moore
#UPDATED: 16 SEP 2022

#AIM: Iterate through values of K, with multiple replicates

#METHODS from MIMS et al 2016:
#Independent putative pops, total of n putative pops
#10 iterations of each K, from 1 to K+1
#500,000 cycles +
#Burn-in of 50,000 cycles +
#Allowed admixture
#Allowed correlated allele frequencies
#LOCPRIOR model (also compared results to without LOCPRIOR)
#detal-K method for most likely K:
##assessed by the second-order rate of change in log-likelihood
##dK cannot be calculated for K=1; K=1 is assumed most likely for runs in which K=1 has greatest log-likelihood, α (Dirichlet paramter for degree of admixture) varies throughout run rather than converging, and assignment to genetic clusters when K>1 tends to be highly admixed w/in indivs
##Identified terminal clusters (K=1) first by looking at log likelihood, second by visual α and individual admixture inspection
#Analysis repeated for clusters where K>1 and n>1 until terminal clusters described.

#notes: before pulling in the proportion of membership table, pull in the lines that read Ln Prob of Data, Mean value of ln likelihood, etc.

#set up run
cd ~/Documents/Projects/console
rm structure_loci.txt
cp ../temporal_aztreefrog/microsat_data/structure_loci.txt ./
sed -i '1d' structure_loci.txt

#make function for loop
iter_struct () {
	local k=$1
	for i in {1..3..1}
	do
		./structure > runK${k}_r$i.txt -m mainparams_testing -K $k -o Year_Pop_K${k}_r$i
        	mv Year_Pop_K${k}_r${i}_f Output_Files/
		mv runK${k}_r$i.txt Run_Files/
	done
	wait
        echo "STRUCTURE with K=$k complete"
}

#run program
#for k in {1..24..1} #seq 1 to 24 (k+1), by 1
for k in 1 2
do
	iter_struct "$k" &
done
wait

echo "STRUCTURE runs complete..." #why isn't it making it to this point??


