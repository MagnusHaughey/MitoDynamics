
#!/bin/bash

inPath=$(readlink -m $1)
inPath=$inPath'/*.tar.gz'

outFileMean=$(readlink -m $1)'/allMean.dat'
outFileVariance=$(readlink -m $1)'/allVariance.dat'


for archive in $(ls -a $inPath)
do


	# Find path to data file in specific format
	filePath=$(tar -ztvf $archive | grep "meanAndVariance.dat" | rev | cut -d" " -f1 | rev)	


	# Re-format mean and variance values and write to single file
	# Mean	
	for dataPoint in $(tar -axf $archive $filePath -O | cut -d" " -f3)
	do
		printf "$dataPoint " >> $outFileMean				
	done

	printf "\n" >> $outFileMean	
	
	# Variance
	for dataPoint in $(tar -axf $archive $filePath -O | cut -d" " -f4)
        do
                printf "$dataPoint " >> $outFileVariance
        done

	printf "\n" >> $outFileVariance
	

done
