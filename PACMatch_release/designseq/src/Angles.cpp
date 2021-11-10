/*
 * Angles.cpp
 *
 *  Created on: 2017Äê11ÔÂ1ÈÕ
 *      Author: notxp
 */

#include "designseq/Angles.h"

namespace NSPdesignseq {

float angle(const XYZ& A, const XYZ& B){
	float x = A*B;
	float cosAng = x/A.length()/B.length();
	if(cosAng > 1)
		cosAng = 1;
	else if(cosAng < -1)
		cosAng = -1;
	return acos(cosAng)*180/3.141592653;
}

float angle(const XYZ& A, const XYZ& B, const XYZ& C)
{
	XYZ BA = A - B;
	XYZ BC = C - B;
	return angle(BA,BC);
}


float dihedral(const XYZ& A, const XYZ& B, const XYZ& C, const XYZ& D)
{
	XYZ BA = A-B;
	XYZ CD = D-C;

	XYZ BC = C-B;
	XYZ CB = B-C;

	XYZ n1 = BA^BC;
	XYZ n2 = CB^CD;

	float ang = angle(n1,n2);

	float symbol = (n1^n2)*BC;
	if(symbol >= 0)
		return ang;
	else
		return 0-ang;
}

} /* namespace NSPdesignseq */
