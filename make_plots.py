
import numpy as np 
import matplotlib.pyplot as plt 
import sys


# Manually put in fit parameters from file 
fit_params = np.loadtxt("/Users/haughe01/Documents/mtDynamics/2020_03_09/cpp/fromApocrita/allVarianceFitParameters.dat", skiprows = 1)


points = []
labels = []
gradients = []
initialFreqs = []

fig = plt.figure(figsize=(11,5))
plt.subplot(1, 2, 1)

colours = ["Blue" , "skyblue" , "teal" , "limegreen" , "darkolivegreen" , "gold" , "darkorange" , "orangered" , "maroon"]

for i in range(1 , len(sys.argv)):

	# Find initial heteroplasmy from filename 
	f_0 = [ word for word in sys.argv[i].split("/") if "originalMutationFrequency=" in word ][0]
	f_0 = f_0.split("=")[-1]

	# Find  fit parameters from file 
	params = [row for row in fit_params if row[0] == float(f_0)][0][1:]

	a = params[0]
	b = params[1]


	initialFreqs.append(float(f_0))
	gradients.append(a)

	generation , variance , err = np.loadtxt(sys.argv[i] , unpack=True)

	lineFit = []
	for j in range(len(generation)):
		lineFit.append(a*j + b)


	scatter = plt.scatter([2**gen for gen in generation] , variance , zorder = 10 , color = colours[i-1] , s = 4 , label = 'Simulation')
	points.append(scatter)
	labels.append(r"$f_0=${}".format(f_0))
	if (float(f_0) >= 0.5):
		line = plt.plot([2**gen for gen in generation] , lineFit , zorder = 5 , c = 'lightblue' , linestyle = ':' , label = '')



plt.xscale("log")

plt.xlabel('Population size', fontsize=20)
plt.ylabel('Var(f)', fontsize=20)


points.append(line[0])
labels.append(r"$f(x) = a*log(x) + b$")
plt.legend([points for points in points] , labels)
plt.xlim([0.5,35000])
plt.ylim(-0.0005,0.008)


#========================================= Plot fit parameters 

a = 0.00194249 

plt.subplot(1, 2, 2)
params_scatter = plt.scatter(initialFreqs , gradients , s = 7)
plt.ylim([0,0.0006])

fitted_line = []

for x in np.linspace(initialFreqs[0]-0.05,initialFreqs[-1]+0.05,1000):
	fitted_line.append(a*x*(1.0-x))

line = plt.plot(np.linspace(initialFreqs[0]-0.05,initialFreqs[-1]+0.05,1000) , fitted_line , color = "lightblue" , zorder = -10 , linestyle='dotted')

plt.legend([params_scatter , line[0]] , ["Simulation" , r"$a(f_0) \propto f_0 \cdot (1 - f_0)$"])

plt.xlabel(r'$f_0$', fontsize=20)
plt.ylabel(r'$a$', fontsize=20)

plt.xticks(np.linspace(0.1,0.9,9))

plt.tight_layout()

#plt.show()

# Export plot
out_fig = 'varianceDataAndFittedParams.png'
plt.savefig(out_fig , dpi=300 , format='png')





















