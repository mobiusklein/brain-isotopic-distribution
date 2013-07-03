/********************************************************************
	created:	2013/07/01
	created:	1:7:2013   0:27
	filename: 	E:\Sync\Dropbox\Codes\BRAIN\src\ConfigLoader.h
	file path:	E:\Sync\Dropbox\Codes\BRAIN\src
	file base:	ConfigLoader
	file ext:	h
	author:		Han Hu
	
	purpose:	Loading config files (xml based).
*********************************************************************/

#ifndef BRAIN_CONFIGLOADER_H
#define BRAIN_CONFIGLOADER_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

namespace brain
{
	class ConfigLoader 
	{
		public:
			virtual void load(const std::string& filename) = 0;
			
			virtual ~ConfigLoader() {}

	};

}
#endif /* BRAIN_CONFIGLOADER_H */
