#include "Peak.h"
#include <iostream>

namespace brain
{
	void Peak::printout(std::ostream& os) const
	{
		std::cout.precision(10);
		os << "MZ: " << mz << "\t"
			<< "Intensity: " << intensity << "\t" << std::endl;
	}
}