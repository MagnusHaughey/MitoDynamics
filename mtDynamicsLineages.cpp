# include <iostream>
# include <fstream>
# include <sstream>
# include <math.h>
# include <dirent.h>
# include <vector>
# include <cstdlib>


using namespace std;




/*******************************************************************************/



// Define poisson distributions
//default_random_engine generator(0);			// give generic seed here, re-seed with user-specified seed later


// Initialise distribution mean to zero, set to non-zero value later
//std::poisson_distribution<int> poisson_t = std::poisson_distribution<int>(0);



/*******************************************************************************/

int mtIndex, numCells, numTagged, numFixated, counter, duplicates, popSize, popSize_previousGen;

const int originalMitoCopyNumber = 1000;
const int numGenerations = 12;



/*******************************************************************************/



// Define a mitochondria 'unit'
class Mitochondria
{

	public:
		vector<int> mutations;


		// Constructor for mitochondria object
		Mitochondria(){}

		// Copy constructor 
		Mitochondria(const Mitochondria &mt)
		{
			mutations = mt.mutations;
		}


		// Set() and get() methods

		void setMutations(vector<int> mutations)
		{
			this->mutations = mutations;
		}

		vector<int> getMutations()
		{
			return this->mutations;
		}


};





// Define a cell
class Cell
{

	public:		
		vector<Mitochondria> allMitochondria;
		//Mitochondria allMitochondria;
		bool tagged;
		int lineage;
		int currentMitoCopyNumber;

		// Constructor for Cell object
		Cell()
		{
			tagged = false;
			lineage = 0;
			currentMitoCopyNumber = 0;
			vector<Mitochondria> allMitochondria(originalMitoCopyNumber , Mitochondria());

			//allMitochondria = temp;
			//allMitochondria.resize(originalMitoCopyNumber);
			//Mitochondria allMitochondria = new Mitochondria[originalMitoCopyNumber];
			//for (int i = 0; i < originalMitoCopyNumber; ++i)
			//{
			//	*allMitochondria[i] = Mitochondria();
			//}
		}

		void setMitochondria(vector<Mitochondria> allMitochondria)
		{
			this->allMitochondria = allMitochondria;
		}

		void addMitochondria( Mitochondria newMito , int index)
		{
			//this->allMitochondria.push_back(newMito);
			this->allMitochondria[index] = newMito;
		}

		void setTag(bool tag)
		{
			this->tagged = tag;
		}

		bool isTagged()
		{
			return this->tagged;
		}

		void setLineage(int lineage)
		{
			this->lineage = lineage;
		}

		int getLineage()
		{
			return this->lineage;
		}
};






int main(int argc, char const *argv[])
{	



	double meanVAF = 0.0;
	double varianceMutatedFreq = 0.0;
	double VAF = 0.0;
	double varSum = 0.0;
	int mutated = 0;
	int startTagging = 6;

	double originalMutationFequency = atof(argv[2]);

	// Create a vector of blank mitochondria to copy into new cells (since I can't find how to do this in the contructor)
	vector<Mitochondria> blankMitochondria(originalMitoCopyNumber , Mitochondria());



	// Seed the random number generator 
	srand48(atoi(argv[1]));



	// Create array to contain each cell in each generation (therefore shape of this array is (numGenerations , 2^{generation})) 
	Cell ** allCells = new Cell*[numGenerations];
	for (int gen = 0; gen < numGenerations + 1; ++gen)
	{

		popSize = (int)pow(2.0 , gen);
		allCells[gen] = new Cell[popSize];
		for (int i = 0; i < popSize; ++i)
		{
			allCells[gen][i] = Cell();
			allCells[gen][i].setMitochondria(blankMitochondria);

			//cout << allCells[gen][i].currentMitoCopyNumber << endl;
		}

	}

	//cout << "Initialised all cells" << endl;





/**************************************************************************************/



	// Create first cell (generation zero)
	Cell firstCell;
	firstCell.setMitochondria(blankMitochondria);



	// Populate first cell with mitochondria, with fraction f = originalMutationFequency mutated
	for (int i = 0; i < originalMitoCopyNumber; ++i)
	{

		//cout << i << endl;
		//firstCell.addMitochondria( Mitochondria() , i );
		//cout << firstCell.allMitochondria.size() << endl;

		if (i < (originalMitoCopyNumber * originalMutationFequency))
		{
			firstCell.allMitochondria[i].mutations.push_back(1);
		}
		else 
		{
			firstCell.allMitochondria[i].mutations.push_back(0);
		}

	}

	firstCell.currentMitoCopyNumber = originalMitoCopyNumber;


	// Add first cell to allCells array (at zeroth generation)
	allCells[0][0] = firstCell;

	//cout << "added first cell" << endl;



	// Create array to hold all information from tagged cells at each generation
	numTagged = (int)(pow( 2 , startTagging ));
	double ** taggingData = new double*[numTagged];

	for (int i = 0; i < numTagged; i++)
	{
		taggingData[i] = new double[numGenerations - startTagging];

		for (int j = 0; j < numGenerations - startTagging; ++j) taggingData[i][j] = 0.0;

	}





/**************************************************************************************/


	// Open data files
	stringstream f;

	f.str("");
	f << "./DATA/numGenerations=" << numGenerations << "_originalMitoCopyNumber=" << originalMitoCopyNumber
			<< "_originalMutationFrequency=" << originalMutationFequency << "/" << atoi(argv[1]);

	DIR *dir = opendir(f.str().c_str());
	if(!dir)
	{
		f.str("");
		f << "mkdir -p ./DATA/numGenerations=" << numGenerations << "_originalMitoCopyNumber=" << originalMitoCopyNumber
			<< "_originalMutationFrequency=" << originalMutationFequency << "/" << atoi(argv[1]);
		system(f.str().c_str());
	}


	ofstream freq_data;
	f.str("");
	f << "./DATA/numGenerations=" << numGenerations << "_originalMitoCopyNumber=" << originalMitoCopyNumber
			<< "_originalMutationFrequency=" << originalMutationFequency << "/" << atoi(argv[1]) << "/meanAndVariance.dat";
	freq_data.open(f.str().c_str());

	ofstream lineage_data;
	f.str("");
	f << "./DATA/numGenerations=" << numGenerations << "_originalMitoCopyNumber=" << originalMitoCopyNumber
			<< "_originalMutationFrequency=" << originalMutationFequency << "/" << atoi(argv[1]) << "/lineages.dat";
	lineage_data.open(f.str().c_str());

	ofstream fixation_data;
	f.str("");
	f << "./DATA/numGenerations=" << numGenerations << "_originalMitoCopyNumber=" << originalMitoCopyNumber
			<< "_originalMutationFrequency=" << originalMutationFequency << "/" << atoi(argv[1]) << "/fixation.dat";
	fixation_data.open(f.str().c_str());


	// Write 'zeroth' generation data to file 
	freq_data << "0 1 " << originalMutationFequency << " 0" << endl;



/**************************************************************************************/



	for (int gen = 1; gen < numGenerations + 1; ++gen)
	{

		popSize = (int)pow(2.0 , gen);			
		popSize_previousGen = pow(2.0 , gen-1);		// population size for previous generation 


		
		// Loop over all cells in 'previous' generation and divide
		for (int i = 0; i < popSize_previousGen; ++i)
		{

			// Parent cell is allCells[gen - 1][i]
			// Daughter cells are allCells[gen][2*i] and allCells[gen][(2*i) + 1]

			
			// Lineage tracking
			if (allCells[gen-1][i].isTagged())	// If parent cell is tagged, pass on tag and lineage number to first daughter cell
			{

				allCells[gen][2*i].setTag(true);
				allCells[gen][2*i].setLineage(allCells[gen-1][i].getLineage());
				allCells[gen-1][i].setTag(false);
			}
			




			// Cell division & mitochondria segregation
			for (int j = 0; j < allCells[gen-1][i].currentMitoCopyNumber; ++j)	// Loop over mitochondria in parent cell
			{


				if (drand48() < 0.5)
				{
					allCells[gen][2*i].addMitochondria( allCells[gen-1][i].allMitochondria[j] , allCells[gen][2*i].currentMitoCopyNumber );
					allCells[gen][2*i].currentMitoCopyNumber += 1;
				}
				else
				{
					allCells[gen][(2*i) + 1].addMitochondria( allCells[gen-1][i].allMitochondria[j] , allCells[gen][(2*i) + 1].currentMitoCopyNumber );
					allCells[gen][(2*i) + 1].currentMitoCopyNumber += 1;
				}

			}




		}



		// Lineage tracking
		if (gen == startTagging) 
		{
			for (int i = 0; i < popSize; ++i) 
			{
				allCells[gen][i].setTag(true);
				allCells[gen][i].setLineage(i);
			}
		}


		//cout << "Got here (gen " << gen << ")" << endl;



		// Replenish mt copy number in all cells
		for (int i = 0; i < popSize; ++i)		// loop over all cells in current generation (i.e. newly created daughter cells)
		{

			//cout << "gen " << gen << " -> population size " << popSize << " checking cell " << i << " (number of mt = " << allCells[gen][i].allMitochondria.size() << ")" << endl;


			while(allCells[gen][i].currentMitoCopyNumber < originalMitoCopyNumber)
			{
				// Randomly select one mt to divide
				mtIndex = (int)(drand48() * (allCells[gen][i].currentMitoCopyNumber-1));

				allCells[gen][i].addMitochondria( allCells[gen][i].allMitochondria[mtIndex] , allCells[gen][i].currentMitoCopyNumber );
				allCells[gen][i].currentMitoCopyNumber += 1;

				/*
				if (allCells[i].allMitochondria[mtIndex].mutations[0] != allCells[i].allMitochondria[allCells[i].allMitochondria.size()-1].mutations[0])
				{
					cout << "Copied MT has different mutations to original MT..........." << endl;
					exit(0);
				}
				*/
				//cout << "Replenising mtDNA copy number (copying mt at index " << mtIndex << ")... -> " << allCells[gen][i].allMitochondria.size() << endl;
			}

		}


		/*
		// Write mutation frequency for each cell in current generation to file
		ofstream gen_file;
		f.str("");
		f << "./DATA/numGenerations=" << numGenerations << "_originalMitoCopyNumber=" << originalMitoCopyNumber
				<< "_originalMutationFrequency=" << originalMutationFequency << "/" << atoi(argv[1]) << "/generation" << gen << ".dat";
		gen_file.open(f.str().c_str());

		for (int i = 0; i < popSize; ++i)
		{

			// Count number of mutated mitochondria for each cell 
			mutated = 0;
			for (int j = 0; j < allCells[gen][i].currentMitoCopyNumber; ++j)
			{
				mutated += allCells[gen][i].allMitochondria[j].mutations[0];
			}

			// Write to file
			gen_file << mutated << " " << allCells[gen][i].currentMitoCopyNumber << endl;


		}

		gen_file.close();
		*/
	

		// After each generation, compute mean and variance of frequency of mutated mitochondria
		// Mean
		mutated = 0;
		for (int i = 0; i < popSize; ++i)
		{
			for (int j = 0; j < allCells[gen][i].currentMitoCopyNumber; ++j)
			{
				mutated += allCells[gen][i].allMitochondria[j].mutations[0];
			}

		}
		//meanVAF = mutated / originalMitoCopyNumber;

		meanVAF = (double)mutated / ((double)popSize * (double)originalMitoCopyNumber);



		// Variance
		varSum = 0.0;
		for (int i = 0; i < popSize; ++i)
		{
			// VAF of mutated mt in individual cell
			mutated = 0;
			VAF = 0.0;
			for (int j = 0; j < allCells[gen][i].currentMitoCopyNumber; ++j)
			{
				mutated += allCells[gen][i].allMitochondria[j].mutations[0];
			}
			VAF = (double)(mutated)/(double)(originalMitoCopyNumber);

			varSum += pow( (VAF - meanVAF) , 2 );
		}

		varianceMutatedFreq = varSum / (double)(popSize - 1);



		// At final generation, check that the same number of cells are tagged as were originally tagged (simply quality check)
		if ((gen == startTagging) || (gen == numGenerations))
		{

			// Count number of cells with a tag
			counter = 0;
			for (int i = 0; i < popSize; ++i) 
			{
				if (allCells[gen][i].isTagged()) counter += 1;
			}


			// Find any cells which have the same tag (this should not happen -- each cell should have a unique tag)
			for (int i = 0; i < numTagged; ++i)
			{
				duplicates = 0;
				for (int j = 0; j < popSize; ++j) 
				{
					if ((allCells[gen][j].isTagged()) && (allCells[gen][j].getLineage() == i)) duplicates += 1;
				}

				if (duplicates > 1)
				{
					cout << "More than one cell has the lineage number: " << i << endl;
				}
			}

			cout << "Generation #" << gen << ": mean = " << meanVAF << ", variance = " << varianceMutatedFreq << ", tagged = " << counter << endl;
		}

		else cout << "Generation #" << gen << ": mean = " << meanVAF << ", variance = " << varianceMutatedFreq << endl;




		// Write data to file
		freq_data << gen << " " << popSize << " " << meanVAF << " " << varianceMutatedFreq << endl;




		// Lineage tracking
		if (gen >= startTagging)
		{
			for (int i = 0; i < numTagged; ++i)	// Loop over total number of cells involved in lineage tracking
			{
				for (int j = 0; j < popSize; ++j) 	// Loop over all cells in system
				{
					if ((allCells[gen][j].isTagged()) && (allCells[gen][j].getLineage() == i))
					{
						mutated = 0;
						for (int k = 0; k < allCells[gen][j].currentMitoCopyNumber; ++k)	// Count number of mutated mitochondria in cell
						{
							mutated += allCells[gen][j].allMitochondria[k].mutations[0];
						}
						taggingData[i][gen - startTagging] = (double)mutated/(double)originalMitoCopyNumber;	// Store the VAF of mutation in array
					}
				}
			}
		}




		// Fixation data
		numFixated = 0;
		for (int i = 0; i < popSize; ++i) 	// Loop over all cells in current generation
		{
			mutated = 0;
			for (int j = 0; j < allCells[gen][i].currentMitoCopyNumber; ++j)	// Count number of mutated mitochondria in cell
			{
				mutated += allCells[gen][i].allMitochondria[j].mutations[0];
			}

			if ((mutated == 0) || (mutated == originalMitoCopyNumber)) numFixated += 1;
		}

		// Write data to file
		fixation_data << gen << " " << numFixated << endl;



		// Delete previous generation 
		//if (gen > 1)
		//{
		//	delete[] allCells[gen-1];
		//}


		//if (gen == 3)
		//{
		//	for (int i = 0; i < popSize; ++i) 	// Loop over all cells in current generation
		//	{
		//		cout << allCells[gen][i].currentMitoCopyNumber << endl;
		//	}

			// Write mt-mutations from first cell to stdout
			//for (int i = 0; i < allCells[gen][popSize - 2].allMitochondria.size(); ++i)
			//{
			//	cout << allCells[gen][popSize - 2].allMitochondria[i].mutations[0] << endl;
			//}
			//exit(0);
		//}



	}



	for (int i = 0; i < numTagged; i++)
	{
		for (int j = 0; j < numGenerations - startTagging; ++j)
		{
			lineage_data << taggingData[i][j] << "\t";
		}
		lineage_data << endl;
	}


	freq_data.close();
	lineage_data.close();
	fixation_data.close();

	

/**************************************************************************************/




	//cout << "After  " << numGenerations << " generations..." << endl;
	//cout << "Number of cells: " << allCells.size() << endl;
	//cout << "Original VAF of mutation: " << originalMutationFequency << endl;
	//cout << "Mean fraction of mutated mitochondria: " << meanVAF << endl; 








	return 0;
}













