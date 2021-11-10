/*
 * geometry_python_module.h
 *
 *  Created on: 2018年5月31日
 *      Author: hyliu
 */

#ifndef GEOMETRY_PYTHON_MODULE_H_
#define GEOMETRY_PYTHON_MODULE_H_
#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/stl_iterator.hpp>
#include "geometry/quatfit.h"

#include <iostream>

template <class T>
inline
std::vector<T> to_std_vector(const boost::python::list & l)
{
    return std::vector<T>(boost::python::stl_input_iterator<T>(l),
                             boost::python::stl_input_iterator< T >( ));
}
// Converts a C++ vector to a python list
template <class T>
inline
boost::python::list to_python_list(std::vector<T> vector) {
    typename std::vector<T>::iterator iter;
    boost::python::list list;
    for (iter = vector.begin(); iter != vector.end(); ++iter) {
        list.append(*iter);
    }
    return list;
}

class QuatFit_python{
public:
	NSPgeometry::QuatFit qfit;
	double rmsd2;
	boost::python::list fitcrd_w(const boost::python::list &ref_coord_list,
			boost::python::list & coord_list,
			boost::python::list weights_list){
		std::vector<double> ref_coord=to_std_vector<double>(ref_coord_list);
		std::vector<double> coord=to_std_vector<double>(coord_list);
		std::vector<double> weights=to_std_vector<double>(weights_list);
		rmsd2=qfit.fitting(ref_coord,coord,weights);
		return to_python_list(coord);
	}
	boost::python::list fitcrd(const boost::python::list &ref_coord_list,
			boost::python::list & coord_list){
		std::vector<double> ref_coord=to_std_vector<double>(ref_coord_list);
		std::vector<double> coord=to_std_vector<double>(coord_list);
		rmsd2=qfit.fitting(ref_coord,coord);
		return to_python_list(coord);
	}

	boost::python::list transform(boost::python::list &coord_list){
		std::vector<NSPgeometry::XYZ> crd;
		NSPgeometry::XYZ::vectortoxyzs(to_std_vector<double>(coord_list),crd);
		qfit.transform(crd);
		std::vector<double> coord;
		NSPgeometry::XYZ::xyzstovector(crd,coord);
		return to_python_list(coord);
	}
};


#endif /* GEOMETRY_PYTHON_MODULE_H_ */
