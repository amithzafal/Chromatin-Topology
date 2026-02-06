#!/bin/bash

#SBATCH --job-name TR
##SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -p Lake
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 02_compute_contact_map_TR.out # File to which STDOUT will be written
#SBATCH -e 02_compute_contact_map_TR.out # File to which STDERR will be written 

cd /scratch/Bio/mdistefano/2024_11_13_Analysis_of_Amith_trajectories/

for simDir in $(ls -1 | grep -v OLD | grep -w TopoNov);
do
    echo $simDir
    cd  $simDir

    for replicaDir in $(ls -1);
    do
	if [[ ! -d $replicaDir ]];
	then
	    continue
	fi
	echo $replicaDir

	for tmax in $(seq 10 10 100);	
	do
	    tmin=$(echo $tmax | awk '{print $1-10}')
	    echo $tmin $tmax

	    check=$(ls -1 /home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/${simDir}/${replicaDir}/time_${tmax}_pos.res 2> /dev/null | wc -l)
	    if [[ $check -lt 1 ]];
	    then
		continue
	    fi
	    #continue
	    
	    pwd
	    bash ~/2024_11_13_Analysis_of_Amith_trajectories/scripts/02_compute_contact_map_PSMN_TR.sh ${simDir} ${replicaDir} ${tmin} ${tmax}
	done # Close cycle over $tmin
    done
    cd .. # Exit simDir
done # Close cycle over $simDir	
wait
