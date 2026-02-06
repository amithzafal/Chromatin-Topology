#!/bin/bash

#SBATCH --job-name TR
##SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -p Lake
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 03_compute_contact_map_TR.out # File to which STDOUT will be written
#SBATCH -e 03_compute_contact_map_TR.out # File to which STDERR will be written 

cd /scratch/Bio/mdistefano/2024_11_13_Analysis_of_Amith_trajectories/
ls -lrtha /scratch/Bio/mdistefano/2024_11_13_Analysis_of_Amith_trajectories/
ls -lrtha /home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/TopoRate_Vary_Extrusion_0/
ls -lrtha /home/aabdulla/3D/SingleChain_July2024/LatticePoly/LatticePoly/data/TopoRate_Vary_Feb2025


for simDir in $1 ; # $(ls -1 | grep -v OLD | grep $1);
do
    simDir=${simDir}_intermingled

    mkdir -p $simDir
    echo $simDir
    cd  $simDir
    #ls -lrtha /home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/TopoRate_Vary_Extrusion_0/${simDir}
    ls -lrtha /home/aabdulla/3D/SingleChain_July2024/LatticePoly/LatticePoly/data/TopoRate_Vary_Feb2025/${simDir}
    pwd

    #for replicaDir in $(ls -1);
    #do
    #if [[ ! -d $replicaDir ]];
    #then
    #    continue
    #fi
    #echo $replicaDir
    
    #for tmax in 100 ; #$(seq 0 1 100);
    for tmax in $(seq 0 1 10);		
    do
	tmin=$(echo $tmax | awk '{print $1-1}')
	echo $tmin $tmax

	inputDir=/home/aabdulla/3D/SingleChain_July2024/LatticePoly/LatticePoly/data/TopoRate_Vary_Feb2025/
	#inputDir=/home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/TopoRate_Vary_Extrusion_0/
	#inputDir=/home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/
	#inputDir=/home/aabdulla/3D/SingleChain_July2024/LatticePoly/LatticePoly/data/
	check=$(ls -1 ${inputDir}/${simDir}/*/time_${tmax}_pos.read_data ${inputDir}/${simDir%_intermingled}/*/time_${tmax}* 2> /dev/null | wc -l)
	echo $check
	if [[ $check -lt 1 ]];
	then
	    ls -1 ${inputDir}/${simDir}/*/time_${tmax}.*
	    continue
	fi
	#continue
	
	pwd
	#bash ~/2024_11_13_Analysis_of_Amith_trajectories/scripts/02_compute_contact_map_PSMN_TR.sh ${simDir} ${replicaDir} ${tmin} ${tmax}
	bash ~/2024_11_13_Analysis_of_Amith_trajectories/scripts/03_compute_contact_map_PSMN_TR.sh ${simDir} all ${tmin} ${tmax} ${inputDir}	
    done # Close cycle over $tmin
    #done # Close cycle over $replicaDir
    cd .. # Exit simDir
done # Close cycle over $simDir	
wait
