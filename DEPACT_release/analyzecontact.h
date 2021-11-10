/*
 * distbins.h
 *
 *  Created on: 2018年8月23日
 *      Author: hyliu
 */

#ifndef ANALYZECONTACT_H_
#define ANALYZECONTACT_H_
#include <vector>
#include <string>
#include <map>
#include <cmath>
namespace NSPgeometry{
class XYZ;
}
namespace myobcode {
void samplerefdist(const std::vector<std::string> &atypes, const std::vector<NSPgeometry::XYZ> &acrd,
		std::map<std::string,std::vector<double>> &distr1,
		std::map<std::string,std::vector<double>> &distr2);
}


#endif /* ANALYZECONTACT_H_ */
