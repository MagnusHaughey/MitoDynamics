



import numpy as np 
import sys

# Import data
inData = np.loadtxt(sys.argv[1])


if ('Lineages' in sys.argv[1]):

	# Generation when tagging begins in original simulation
	startTagging = 6

	outPath = sys.argv[1] + ".processed.dat"
	f = open(outPath , 'w')

	# Compute mean and variance for each generation
	for i in range(inData.shape[1]):
		f.write("{} {} {}\n" .format( i+startTagging , np.mean(inData[:,i]) , np.var(inData[:,i]) ))


else:

	outPath = sys.argv[1] + ".averaged.dat"
	f = open(outPath , 'w')

	for i in range(inData.shape[1]):
		f.write("{} {} {}\n" .format( i , np.mean(inData[:,i]) , np.std(inData[:,i])/np.sqrt(inData.shape[0]) ))

f.close()