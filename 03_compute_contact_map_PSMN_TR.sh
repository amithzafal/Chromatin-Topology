# Arguments

echo "wDir: $PWD"

simDir=$1

modelResolution=1
mapResolution=1

res=1 # 1 bead
ncopies=1 # In this way the 2 chains won't be considered as two copies of the same chromosome

nCopies=3
b=1.0
dir=$(echo ${PWD} | sed "s,/, ,g" | awk '{print $NF}')

for rc in 1.5 ; # in lattice-units, the same units of the coordinates file
do
    echo "rc (nm) = ${rc}"
    for tmin in ${3} ;
    do
	for tmax in ${4} ;
	do
	    if [[ ${tmin} -ge ${tmax} ]];
	    then
		echo "t-min ${tmin} is larger than t-max ${tmax}"
		continue
	    fi	    
	    for tdelta in 1 ;		      
	    do
		echo "timestep: From ${tmin} to ${tmax} by ${tdelta}"	
		
		for copies in intra ;
		do
		    dirName=$(echo $dir | sed -e "s/\./_/g")
		    echo "dirname : $dirName"
		    
		    #for replicaDir in ${2};
		    #do
		    nparticles=$(echo 2000 $nCopies | awk '{print $1*$2}') # All the particles in the system
		    echo ${res} ${nparticles} ${rc} ${b} ${nCopies}		    
		    #if [[ ! -d ${replicaDir} ]];
		    #then
		    #echo "${replicaDir} doesn't exist!"
		    #continue
		    #fi
		    
		    #outMatrix=rc_${rc}nm_res_${mapResolution}bead_from_${tmin}_to_${tmax}_every_${tdelta}_${replicaDir}.tab
		    outMatrix=rc_${rc}nm_res_${mapResolution}bead_at_${tmax}_${simDir}.tab
		    i=0
		    echo "outMatrix : ${outMatrix}"
		    
		    tmpDir=_tmp_rc_${tmax}_TR
		    if [[ -d ${tmpDir} ]];
		    then
			echo "${tmpDir} exists!"
			continue
		    fi
		    
		    outFile=${outMatrix%.tab}_TR.tab
		    ls -lrtha ${outFile} 2> /dev/null
		    
		    if [[ ! -e ${outFile} ]];
		    then
			mkdir ${tmpDir}
			
			cd ${tmpDir}
			touch ${outMatrix}
			
			rm -fvr DATA_FILES_INPUT.txt			
			echo ${tmin} ${tdelta} ${tmax}
			for t in $tmax ; #$(seq ${tmin} ${tdelta} ${tmax});
			do
			    echo $t
			    #for r in replicaDir ;
			    #do
			    inputDir=$5
			    #inputDir=/home/aabdulla/3D/March2024/LatticePoly/LatticePoly/data/
			    #inputDir=/home/aabdulla/3D/SingleChain_July2024/LatticePoly/LatticePoly/data/

			    ls -1 ${inputDir}/${simDir%_intermingled}/*/*_${t}_pos.re* 2> /dev/null > DATA_FILES_INPUT.txt

			    #done # Close cycle over ${r}
			done # Close cycle over ${t}		       
			cat DATA_FILES_INPUT.txt
			
			nlines=$(cat DATA_FILES_INPUT.txt | wc -l | awk '{print $1}')
			if [[ $nlines -lt 11 ]];
			then
			    rm -fr ${outMatrix} ${outMatrix%.tab}.png output.log distances.txt DATA_FILES_INPUT.txt	    
			    continue
			fi
			echo "Computing the model matrix using ${nlines} snapshots"
			
			if [[ ! -e _chromosomes ]];
                        then
                            rm -fvr _chromosomes
                            for nc in $(seq 1 1 ${nCopies});
                            do
				offset=$(echo $nc 2000 | awk '{print int((($1-1)*$2))}')
                                awk -v nc=$nc 'BEGIN{for(i=0;i<=40;i++){chr[i]=i}; s1=(nc-1)*2000; s2=nc*2000; for(i=s1;i<s2;i++){print i,chr[nc-1]}}' >> _chromosomes
                            done    
			fi
			#cat _chromosomes
			
			echo $outMatrix #$outProb
			touch ${outMatrix}
			#touch ${outFile}
			
			# -r:		 Resolution of the map in monomers
			# -p:		 Number of particles in the system
			# -d:		 Contact distance cutoff between particles (nm)
			# -b:		 Monomer diameter (nm)
			# -k:		 Ploidy of the system
			# -i:                Consider (0) or not (1) inter copies contacts
			echo ${res} ${nparticles} ${rc} ${b} ${ncopies} 0
			np=$(echo ${nparticles} | awk '{print $1}')
			~/2024_11_13_Analysis_of_Amith_trajectories/scripts/compute_contact_map_TR -r ${res} -p ${np} -d ${rc} -b ${b} -k ${ncopies} -i 0
			#i, chrom1, start1, end1, compName, (int) cis[i], Ccis[i], (int) trans[i], Ctrans[i], (int) tot[i], Ctot[i], cis[i], trans[i]);
			mv TR_per_bin.txt ${outFile}
			
			mv contacts.tab ${outMatrix}
			rm -fr output.log distances.txt DATA_FILES_INPUT.txt contacts.tab
			
			#rm ${outMatrix}
			rsync -avz ${outMatrix%.tab}_TR.tab ../
			cd ..
		    fi
		    
		    # Compute TR per chromosome
		    inFile=${outFile}
		    outFile=${outMatrix%.tab}_TR_perChrom.tab
		    if [[ ! -e ${outFile} ]];
		    then
			cd ${tmpDir}
			rsync -avz ../${inFile} .
			
			# TR all
			echo "#chrom cisChromContactsFraction transChromContactsFraction" | awk '{for(i=1;i<NF;i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFile}
			for chrom in $(seq 0 1 $((${nCopies}-1)));
			do
			    echo $chrom
			    sed "s/_/ /g" ${inFile} | grep -w $chrom | awk -v c=${chrom} '{if($2==c){cis+=$(NF-1); trans+=$NF; cnt++}}END{if(cnt!=0){print c,cis/cnt,trans/cnt}}' >> ${outFile}
			    
			    echo ${1}
			    cat ${outFile}
			done # Close cycle over ${chrom}
			rsync -avz ${outMatrix%.tab}_TR_perChrom.tab ../
			cd ..
		    fi
		    rm -fvr ${tmpDir}
		    #exit
		    #done # Close cycle over ${replicaDir}
		done # Close cycle over ${copies}			
	    done # Close cycle over ${nreplicas}		    		
	done # Close cycle over ${tmin}
    done # Close cycle over ${tmax}
done # Close cycle over ${rc}
