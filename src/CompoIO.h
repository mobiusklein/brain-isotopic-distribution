/********************************************************************
	created:	2013/07/02
	created:	2:7:2013   12:35
	filename: 	CompoIO.h
	file path:	BRAIN\src
	file base:	CompoIO
	file ext:	h
	author:		Han Hu (hh.earlydays@gmail.com)
	
	purpose:	Input and output based on specified format.
*********************************************************************/

#ifndef BRAIN_COMPOIO_H
#define BRAIN_COMPOIO_H

#include "Composition.h"
#include "PeakList.h"
#include "Param.h"

namespace compo_io
{
	using namespace std;
	using namespace param;
	using namespace brain;
	class CompoIO
	{
	public:
		// Constructor.
		CompoIO()
			: param(Param::Instance())
		{}

		// Input format:
		vector<pair<Composition, int> > readCSVFile(const string& filename);

		// Output format:
		void exportDistribution(ostream& os, const Composition& compo, AggregatedIsotopicVariants peaklist);

	private:
		//std::vector<Composition> compo_vec;
		Param& param;
	};
}

#endif /* BRAIN_COMPOIO_H */