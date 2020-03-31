

total=0
seed=$1
initialHeteroplasmy=$2
dseed=$3

while :
do
	
	# Simulate data
	printf "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
	./mtDynamicsLineages $seed $initialHeteroplasmy

	
	# Compress data
	filePath='DATA/numGenerations=16_originalMitoCopyNumber=1000_originalMutationFrequency='"$initialHeteroplasmy"'/'"$seed"
	tar czf $(echo $filePath'.tar.gz') $filePath
	rm -r $filePath


	seed=$(( seed + dseed ))
	total=$(( total + 1 ))


	if (( total == 20 ))
	then

		break

	fi

done
