/*
 * Angles.h
 *
 *  Created on: 2017Äê11ÔÂ1ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ANGLES_H_
#define DESIGNSEQ_ANGLES_H_
#include "geometry/xyz.h"

namespace NSPdesignseq {

using namespace NSPgeometry;

float angle(const XYZ& A, const XYZ& B);
float angle(const XYZ& A, const XYZ& B, const XYZ& C);
float dihedral(const XYZ& A, const XYZ& B, const XYZ& C, const XYZ& D);


} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ANGLES_H_ */
