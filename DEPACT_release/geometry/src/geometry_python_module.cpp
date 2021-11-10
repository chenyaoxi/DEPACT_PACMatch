/*
 * geometry_python_module.cpp
 *
 *  Created on: 2018年5月31日
 *      Author: hyliu
 */
#include "geometry/geometry_python_module.h"
BOOST_PYTHON_MODULE(molgeometry)
{
	using namespace boost::python;
	using namespace NSPgeometry;
	class_<QuatFit_python>("QuatFit",init<>())
			.def("fitcrd",&QuatFit_python::fitcrd)
			.def("fitcrd_w",&QuatFit_python::fitcrd_w)
			.def("transform",&QuatFit_python::transform)
			.def_readonly("rmsd2",&QuatFit_python::rmsd2);
}






