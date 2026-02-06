#!/bin/bash
1;95;0c
#SBATCH --job-name cisCounts
##SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -p Lake
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o cisCounts.out # File to which STDOUT will be written
#SBATCH -e cisCounts.out # File to which STDERR will be written

# Prepare micro-C datasets
echo "time totCounts cisCounts transCounts" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > cisCounts_exp.tab


mainDir=/scratch/Bio/mdistefano/2024_11_13_Analysis_of_Amith_trajectories/
mkdir -p ${mainDir}/cisCounts_distributions
mkdir -p ${mainDir}/transCounts_distributions

tMax=20

#for condition in ExtNov ExtTopoNov NoneNov TopoNov ;
#for condition in TopoNov ;
#for condition in Ext1000_Feb2025 None_Feb2025 TopoExt_Feb2025 Topo_Feb2025 ;
for condition in  TopoRate_1.0000e-01 TopoRate_1.0000e-01_intermingled TopoRate_1.0000e-02 TopoRate_1.0000e-02_intermingled TopoRate_1.0000e-03 TopoRate_1.0000e-03_intermingled TopoRate_1.0000e-04 TopoRate_1.0000e-04_intermingled TopoRate_1.0000e-05 TopoRate_1.0000e-05_intermingled TopoRate_1.0000e-06 TopoRate_1.0000e-06_intermingled TopoRate_1.0000e-07 TopoRate_1.0000e-07_intermingled ;
do
    if [[ ! -d ${mainDir}/${condition} ]];
    then
	ls -lrtha ${mainDir}/${condition}
	continue
    fi
    echo $condition
    outFile=${mainDir}/cisCounts_distributions/cisCounts_over_time_tMax_${tMax}_sameScale_${condition}.pdf

    
    if [[ ! -e ${outFile%pdf}tab ]];
    then
	cp cisCounts_exp.tab cisCounts.tab
	
	#for t in $(seq 0 1 100) ;
	for t in $(seq 0 1 ${tMax}) ;		 
	do
	    th=$(echo $t) ; echo $th
	    for rc in 1.5 ;
	    do
		#cat ${mainDir}/${condition}/rc_${rc}nm_res_1bead_from_*_to_${t}_every_1_*_TR.tab  | grep -v "#" | awk -v d=$condition -v rc=${rc} -v t=${th} '{print d"_"t,$(NF-2),$3,$5}'   | sed -e "s/Nov//g" >> cisCounts.tab
		cat ${mainDir}/${condition}/rc_${rc}nm_res_1bead_at_${t}_*_TR.tab  | grep -v "#" | awk -v d=$condition -v rc=${rc} -v t=${th} '{print t,$(NF-2),$3,$5}'   | sed -e "s/Nov//g" >> cisCounts.tab
	    done # Close cycle over $rc
	done # Close cycle over $t

	awk '{print $1}' cisCounts.tab | uniq
	wc -l cisCounts.tab
	cp cisCounts.tab ${outFile%pdf}tab
    fi
    cp ${outFile%pdf}tab cisCounts.tab 
    
    # Make the plot
    if [[ ! -e ${outFile} ]];
    then
	sed -e "s/#${condition}//g" scripts/04_make_cisCounts_violinPlots.R > _tmp.R
	Rscript _tmp.R 1860
	rm -fvr _tmp.R
	
	mv -v cisCounts.pdf ${outFile}
    fi

    if [[ ! -e ${mainDir}/transCounts_distributions/transCounts_over_time_tMax_${tMax}_sameScale_${condition}.pdf ]];
    then
	sed -e "s/#${condition}//g" scripts/04_make_transCounts_violinPlots.R > _tmp.R
	Rscript _tmp.R 1800
	rm -fvr _tmp.R
	
	mv -v transCounts.pdf ${mainDir}/transCounts_distributions/transCounts_over_time_tMax_${tMax}_sameScale_${condition}.pdf
    fi
    rm cisCounts.tab

done # Close cycle over $condition
