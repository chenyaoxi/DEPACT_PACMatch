/*
 * move3d.h
 *
 *  Created on: 2018年8月8日
 *      Author: hyliu
 */

#ifndef MOVE3D_H_
#define MOVE3D_H_
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
namespace subsitedesign {
/**
 * @brief a rigid body transformation
 *
 * @TODO source file should be moved to another location
 */
struct Move3D{
	double rmatrix[3][3];
	std::vector<double> tvec{0.0,0.0,0.0};
	template<class Archive>
	    void serialize(Archive & ar, const unsigned int version)
	    {
			for(int i=0;i<3;i++)
				for(int j=0;j<3;j++) ar & rmatrix[i][j];
			ar & tvec;
	    }
	std::vector<double> move(const std::vector<double> &point) const {
		std::vector<double> result(3,0.0);
		for(int i=0;i<3;++i){
			double &ri=result[i];
			for(int j=0;j<3;++j){
				ri +=rmatrix[i][j]*point[j];
			}
		}
		for(int i=0;i<3;++i) result[i]+=tvec[i];
		return result;
	}
	std::string to_string() const {
		std::stringstream oss;
		for(int i=0;i<3;++i){
			for(int j=0;j<3;++j) oss <<std::setw(11)<<std::setprecision<<rmatrix[i][j];
			oss<<std::endl;
		}
		for(int j=0;j<3;++j) oss <<std::setw(11)<<std::setprecision<<tvec[j];
		oss<<std::endl;
		return oss.str();
	}
};
}


#endif /* MOVE3D_H_ */
