#include "Param.h"
#include <iostream>


namespace param
{
	void Param::load( const std::string& filename )
	{
		std::cout << "Load " << filename << std::endl;

		using boost::property_tree::ptree;
		ptree pt;

		read_xml(filename, pt);

		BOOST_FOREACH(ptree::value_type &v, pt.get_child("parameters.SystemParameters"))
		{
			if(v.first == "Parameter")
			{
				this->setParameter(v.second.get<std::string>("Name"), 
					v.second.get<std::string>("Value"));
			}
		}

		pt.erase("parameters.SystemParameters");
	}

}