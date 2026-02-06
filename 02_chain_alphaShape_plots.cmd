#!/bin/bash

#SBATCH --job-name alphaShape
#SBATCH -n 1                   # Number of cores. For now 56 is the number max of core available
#SBATCH -p Lake
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o chain_alphaShape_plots.out # File to which STDOUT will be written
#SBATCH -e chain_alphaShape_plots.out # File to which STDERR will be written

cd /scratch/Bio/mdistefano/2024_11_13_Analysis_of_Amith_trajectories/
ls -lrtha

#for dir in $(ls -1 | grep -i "july\|nov" | grep -v OLD);
#for dir in $(ls -1 | grep -wi "nonenov" | grep -v OLD);
#for dir in $(ls -1 | grep -wi "toponov" | grep -v OLD);
#for dir in $(ls -1 | grep -wi "exttoponov" | grep -v OLD);
for dir in $(ls -1 | grep -wi $1 | grep -v OLD);
do
    cd $dir
    pwd
    ls -1 | tail

    if [[ ! -e chain_alphaShape_volume.tab ]];
    then
	rm -fvr chain_alphaShape_volume.tab chain_particlesInHeteroShape.tab
	for time in $(seq 0 1 100);
	do
	    for chain in A B C ;
	    do

		#ls -lrtha ./Ext1000A0135/*/?volume.tab
		#head ./Ext1000A0135/*/?volume.tab
		#ls -lrtha ./Ext1000A0135/*/?particles_in_*.tab | wc -l
		#head ./Ext1000A0135/*/?particles_in_*.tab | wc -l
		#awk -v c=${chain} -v t=${time} '{printf("%s%s\t%d\t",$1,c,t); for(i=2;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' ./*/time_${time}_pos/${chain}volume.tab >>  chain_alphaShape_volume.tab
		awk -v t=${time} '{printf("%s\t%d\t",$1,t); for(i=2;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' ./*/time_${time}_pos*/${chain}volume.tab >>  chain_alphaShape_volume.tab
		tail  chain_alphaShape_volume.tab
		wc -l chain_alphaShape_volume.tab
		
		#awk -v c=${chain} -v t=${time} '{printf("%s%s\t%d\t",$1,c,t); for(i=2;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' ./*/time_${time}_pos/${chain}particles_in_*.tab >>  chain_particlesInHeteroShape.tab
		awk -v t=${time} '{printf("%s\t%d\t",$1,t); for(i=2;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' ./*/time_${time}_pos*/${chain}particles_in_*.tab >>  chain_particlesInHeteroShape.tab
		tail  chain_particlesInHeteroShape.tab
		wc -l chain_particlesInHeteroShape.tab
		
	    done # Close cycle over $chain
	done # Close cycle over $time
    fi
    Rscript ~/2024_11_13_Analysis_of_Amith_trajectories/scripts/02_chain_alphaShape_plots.R	    

    cd ..
done # Close cycle over $dir



