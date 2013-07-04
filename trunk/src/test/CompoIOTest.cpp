/********************************************************************
	created:	2013/07/02
	created:	2:7:2013   11:17
	filename: 	IsotopicDistributionTest.cpp
	file path:	BRAIN\src\test
	file base:	IsotopicDistributionTest
	file ext:	cpp
	author:		Han Hu
	
	purpose:	The test file for brain algorithm.  
						Usage: brain.exe filename
*********************************************************************/

#include "IsotopicDistribution.h"
#include "CompoIO.h"
#include <time.h>
#include <boost/filesystem.hpp>

#pragma warning(disable : 4996)

// argv[1] -- input file.
using namespace std;
using namespace brain;
using namespace compo_io;
using namespace boost::filesystem;

int main(int argc, char* argv[])
{
	try {
		string filebase;
		if(argv[1] == NULL) {
			cout << "No input file, use example input instead!" << endl;
			filebase = "./data/test_input.csv";
		} else
			filebase = argv[1];

		CompoIO compo_io;
		string infile(filebase);
		
		vector<pair<Composition, int> > compo_vec = compo_io.readCSVFile(infile);

		clock_t t;

		//string outfile(filebase);
		// Parse the input path and generate output file name.
		path p(filebase);
		path outfile = p.parent_path();
		string stem = p.stem().string() + "_result";

		outfile /= stem.append(p.extension().string());

		//outfile.append("_result.csv");
		std::filebuf fb;
		fb.open(outfile.c_str(), std::ios::out);
		ostream os(&fb);

		for(vector<pair<Composition, int> >::iterator iter = compo_vec.begin(); iter != compo_vec.end(); iter++)
		{
			// Initial time.
			t = clock();

			IsotopicDistribution iso_dist(iter->first, iter->second);
			AggregatedIsotopicVariants peakset = iso_dist.getAggregatedIsotopicVariants();
			double avg_mass = iso_dist.getAverageMass();

			// Ending time.
			t = clock()-t;
			cout << "It took " << ((double)t)/CLOCKS_PER_SEC << " seconds to calculate " << iter->first.getCompositionString() << endl;

			compo_io.exportDistribution(os, iter->first, peakset, avg_mass);
			
		}
		fb.close();
	}
	catch(std::exception& e)
	{
		cout << "Exception: " << e.what() << endl;
	}

	cout << "The End!" << endl;
	//cin.get();

	return 0;
}