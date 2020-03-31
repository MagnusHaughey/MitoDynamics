
#!/bin/bash

inPath=$(readlink -m $1)
inPath=$inPath'/*.tar.gz'

outPath=$(readlink -m $1)'/allLineages.dat'

for archive in $(ls -a $inPath)
do

	# Find path to data file in specific format
	filePath=$(tar -ztvf $archive | grep "lineages.dat" | rev | cut -d" " -f1 | rev)	

	tar -axf $archive $filePath -O >> $outPath 

done
