#!/bin/bash

#SBATCH --job-name aShape3D
#SBATCH -n 1                   # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -p Lake
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o aShape3D.out # File to which STDOUT will be written
#SBATCH -e aShape3D.out # File to which STDERR will be written

#/home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/TopoNov/
#/home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/NoneNov/
#/home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/ExtNov/
#/home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/ExtTopoNov/

#inputDir=/home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/
#inputDir=/home/aabdulla/3D/SingleChain_July2024/LatticePoly/LatticePoly/data/
inputDir=/home/aabdulla/3D/SingleChain_July2024/LatticePoly/LatticePoly/data/TopoRate_Vary_Feb2025/
#inputDir=/home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/TopoRate_Vary_Extrusion_0/


workingDir=/scratch/Bio/mdistefano/2024_11_13_Analysis_of_Amith_trajectories/
mkdir -p ${workingDir}

#for wDir in TopoJuly NoneJuly ExtJuly ;
#for wDir in TopoNov NoneNov ExtNov ExtTopoNov ;
#for wDir in Topo_Feb2025 None_Feb2025 Ext1000_Feb2025 TopoExt_Feb2025
for wDir in $1"_intermingled"
do
    mkdir -p ${workingDir}/${wDir}/
    cd ${workingDir}/${wDir}/

    for inDir in $(ls -1 ${inputDir}/${wDir%_intermingled} | grep -v res | grep rep | grep -v gz | grep -v tar);
    do
	echo $inDir
	if [[ -d $inDir ]];
	then
	    echo "${inDir} already analyzed"
	    continue
	fi

	mkdir -p $inDir
	cd $inDir
	
	#ls -lrtha ${inputDir}/${wDir}/$inDir
	#for inFile in $(ls -1 ${inputDir}/${wDir}/${inDir}/time_*_pos.res);
	for inFile in $(ls -1 ${inputDir}/${wDir%_intermingled}/${inDir}/time_*_data);	
	do
	    echo $inFile

	    inFile=$(echo ${inFile} | sed -e "s/\// /g" | awk '{print $NF}')
	    outDir=$(echo ${inFile} | sed -e "s/\.res//g" -e "s/\// /g" | awk '{print $NF}')
	    echo ${outDir}
	    if [[ -e ${outDir}/Avolume.tab ]];
	    then
		continue
	    fi
	    
	    echo $outDir
	    mkdir -p ${outDir}
	    cd ${outDir}
	    pwd
	    echo $inFile
	    cat ${inputDir}/${wDir%_intermingled}/${inDir}/${inFile} | awk '{if(NR>=1    && NR<=2000){print $1,$2,$3}}' > _xyzAall
	    cat ${inputDir}/${wDir%_intermingled}/${inDir}/${inFile} | awk '{if(NR>=2001 && NR<=4000){print $1,$2,$3}}' > _xyzBall	     
	    cat ${inputDir}/${wDir%_intermingled}/${inDir}/${inFile} | awk '{if(NR>=4001 && NR<=6000){print $1,$2,$3}}' > _xyzCall
	    wc -l _xyz?all

	    Rscript ~/2024_11_13_Analysis_of_Amith_trajectories/scripts/01_compute_ashape3D_volumeNparticles.R
	    cd .. # Exit outDir
	done
	cd .. # Exit inDir
	#exit	
    done
done
exit
