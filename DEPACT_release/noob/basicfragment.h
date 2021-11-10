/*
 * basicfragment.h
 *
 *  Created on: 2018年8月9日
 *      Author: hyliu
 */

#ifndef BASICFRAGMENT_H_
#define BASICFRAGMENT_H_
#include <string>
#include <vector>
#include <map>
#include <bitset>
#include <dataio/datapaths.h>
using namespace NSPdataio;
namespace subsitedesign {
/**
 * @brief similar to a functional group. Used to define minimum substructures in a ligand
 *
 * It is a user-defined SMARTS pattern read from an input file
 */
struct BasicFragment {
	/**
	 * @brief return the user-defined basic fragments read from an input
	 */
	static const std::map<std::string, BasicFragment> & getmap(const std::string &
			filename = getenvpath("DEPACT_DATAPATH")+"basicfragments.txt") {
		static std::map<std::string, BasicFragment> map;
		static bool initialized{false};
		if(!initialized){
			initmap(&map,filename);
			initialized=true;
		}
		return map;
	}
	std::string smarts;   ///<SMARTS string
	static void initmap(std::map<std::string, BasicFragment> *fragmentmap,
			const std::string &filename);
};
}

#endif /* BASICFRAGMENT_H_ */
