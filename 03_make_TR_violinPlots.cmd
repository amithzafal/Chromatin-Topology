#!/bin/bash

#SBATCH --job-name TR
#SBATCH -n 1                   # Number of cores. For now 56 is the number max of core available
#SBATCH -p Lake
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 03_make_TR_violinPlots.out # File to which STDOUT will be written
#SBATCH -e 03_make_TR_violinPlots.out # File to which STDERR will be written

# Prepare micro-C datasets
echo "sample chrom avg_cis_ratio avg_trans_ratio" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > AvgTransRatio_perChrom_exp.tab


mainDir=/scratch/Bio/mdistefano/2024_11_13_Analysis_of_Amith_trajectories/
mkdir -p ${mainDir}/TR_distributions

#for condition in ExtNov ExtTopoNov NoneNov TopoNov ;
for condition in Ext1000_Feb2025 None_Feb2025 TopoExt_Feb2025 Topo_Feb2025 ;
#for condition in $1 ;
do
    if [[ ! -d ${mainDir}/${condition} ]];
    then
	continue
    fi
    echo $condition
    outFile=${mainDir}/TR_distributions/AvgTransRatio_perChrom_${condition}.pdf
    if [[ -e ${outFile} ]];
    then
	continue
    fi
    
    cp AvgTransRatio_perChrom_exp.tab AvgTransRatio_perChrom.tab
    #awk '{print $1}' AvgTransRatio_perChrom.tab | uniq   
    #rm -fvr AvgTransRatio_perChrom.tab
    
    for t in $(seq 10 10 100) ;    
    do
	th=$(echo $t) ; echo $th
	for rc in 1.5 ;
	do
	    #rc_1.5nm_res_1bead_at_97_None_Feb2025_TR_perChrom.tab
	    cat ${mainDir}/${condition}/rc_${rc}nm_res_1bead_from_*_to_${t}_every_1_*_TR_perChrom.tab  | grep -v "#" | awk -v d=$condition -v rc=${rc} -v t=${th} '{print d"_"t,$(NF-2),$(NF-1),$NF}'   | sed -e "s/Nov//g" >> AvgTransRatio_perChrom.tab
	    #tail AvgTransRatio_perChrom.tab
	done # Close cycle over $rc
    done # Close cycle over $t
    
    awk '{print $1}' AvgTransRatio_perChrom.tab | uniq
    wc -l AvgTransRatio_perChrom.tab

    # Make the plot
    sed -e "s/#${condition}//g" scripts/03_make_TR_violinPlots.R > _tmp.R
    Rscript _tmp.R
    rm -fvr _tmp.R
    
    mv -v AvgTransRatio_perChrom.pdf ${outFile}
    mv -v AvgTransRatio_perChrom.tab ${outFile%pdf}tab
    #exit
done # Close cycle over $condition
