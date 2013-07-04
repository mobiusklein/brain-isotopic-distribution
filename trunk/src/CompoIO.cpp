#include "CompoIO.h"
#include <sstream>
#include <fstream>
#include <boost/xpressive/xpressive.hpp>
#include <boost/lexical_cast.hpp>

namespace compo_io
{
	
	vector<pair<Composition, int> > CompoIO::readCSVFile( const string& filename )
	{
		vector<pair<Composition, int> > compo_vec;

		// Read composition list from csv file.
		ifstream compo_file;
		
		compo_file.open(filename.c_str());

		using namespace boost::xpressive;

		if(compo_file.is_open()) {
			while(!compo_file.eof()) {
				string ROW;
				// Read file line by line.
				getline(compo_file, ROW);

				// Construct regular expression.
				sregex res = bos >> (s1 = +(+_w >> *_d)) >> +_s >> (s2 = *_d);

				// Iterate over all non-white space in the row string.
				//sregex_token_iterator begin(ROW.begin(), ROW.end(), res, -1), end;
				smatch what;
				if(regex_match(ROW, what, res)) {
					
					std::string iso_variants = what[2];

					int num_iso_variants = iso_variants.empty() ? param.getParameter<int>("num_iso_variants").first : boost::lexical_cast<int>(what[2]);

					// Store it into the member variable if molecular composition.
					compo_vec.push_back(make_pair(Composition(what[1]), num_iso_variants));
					
				}

			}

			compo_file.close();
		} else {
			throw runtime_error("Unable to open file!");
		}

		return compo_vec;
	}

	void CompoIO::exportDistribution( ostream& os, const Composition& compo, AggregatedIsotopicVariants peaklist, double avg_mass)
	{
		// First line: Composition
		os << compo.getCompositionString() << "\n";
		
		// Second line: Monoisotopic mass.
		os << setprecision(15);
		os << compo.getMass() << "\n";
		
		// Third line: Average mass.
		os << compo.getAverageMass() << "\n";

		// Fourth line: Average mass based on the value returned by isotopic peaks.
		os << avg_mass << "\n";

		// Theoretical isotopic distribution.
		
		peaklist.printPeakList<peak_mz>(os);

		os << "\n";

	}



}