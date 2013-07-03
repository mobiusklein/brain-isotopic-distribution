#include "Element.h"

namespace brain
{
	
	double Element::getAverageMass() const
	{
		double avg_mass = 0.0;
		const IsotopeSetSequential& iso_index = isotopes.get<0>();

		for(IsotopeSetSequential::const_iterator iter = iso_index.begin(); iter != iso_index.end(); iter++)
		{
			avg_mass += iter->mass * iter->abundance;
		}

		return avg_mass;
	}

}