#include "CompoIO.h"
#include <boost/make_shared.hpp>

// argv[0] -- The filename for input.
// argv[1] -- The filename for output.

int main(int argc, char* argv[])
{
	using namespace std;
	using namespace brain;
	using namespace compo_io;

	CompoIO compo_reader;

	string infile = "./data/test_input.csv";
	string outfile = "./data/test_output.csv";

	vector<pair<Composition, int> > compo_vec = compo_reader.readCSVFile(infile);

	AggregatedIsotopicVariants peaklist;

	PeakPtr pk9 = boost::make_shared<Peak>(471.2, 1320.5);
	peaklist.addPeak(pk9);

	PeakPtr pk10 = boost::make_shared<Peak>(472.2, 170.6);
	peaklist.addPeak(pk10);

	PeakPtr pk11 = boost::make_shared<Peak>(473.2, 222222.7);
	peaklist.addPeak(pk11);

	PeakPtr pk12 = boost::make_shared<Peak>(474.2, 26.4);
	peaklist.addPeak(pk12);

	PeakPtr pk13 = boost::make_shared<Peak>(475.2, 32.5);
	peaklist.addPeak(pk13);

	std::filebuf fb;
	fb.open(outfile.c_str(), std::ios::out);
	ostream os(&fb);
	compo_reader.exportDistribution(os, compo_vec.at(0).first, peaklist);

	fb.close();

	cout << "The End!" << endl;
	cin.get();

	return 0;
}