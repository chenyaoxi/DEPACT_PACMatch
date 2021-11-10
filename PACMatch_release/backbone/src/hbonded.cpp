/*
 * hbonded.cpp
 *
 *  Created on: 2016年12月10日
 *      Author: hyliu
 */

#include "backbone/hbonded.h"

using namespace NSPproteinrep;
using namespace NSPgeometry;

unsigned int NSPproteinrep::hbonded(const BackBoneSite & bs1,
		const BackBoneSite &bs2, HBondGeometry *hbg) {
	if (n1o2hbonded(bs1, bs2, hbg))
		return 1;
	else if (n1o2hbonded(bs2, bs1, hbg))
		return 2;
	return 0;
}

bool NSPproteinrep::n1o2hbonded(const BackBoneSite & bs1,
		const BackBoneSite &bs2, HBondGeometry *hbg) {
	const double RHAMIN = 1.5;
	const double RHAMAX = 2.45;
	const double RDAMAX = 3.5;
	const double RDAMIN = 2.3;
	const double TDHAMIN = 120.0;
	const double THABMIN = 110.0;
	const double TDABMIN = 110.0;
	if (bs1.resname == "PRO")
		return false;
	XYZ n1crd = bs1.getcrd(BackBoneSite::NCRD);
	XYZ h1crd = bs1.hcrd();
	XYZ o2crd = bs2.getcrd(BackBoneSite::OCRD);
	double dh1o2 = (h1crd - o2crd).squarednorm();
	double dn1o2 = (n1crd - o2crd).squarednorm();
	if (dh1o2 > dn1o2)
		return false;
	if (!(dh1o2 > RHAMIN * RHAMIN && dh1o2 < RHAMAX * RHAMAX))
		return false;
	if (!(dn1o2 > RDAMIN * RDAMIN && dn1o2 < RDAMAX * RDAMAX))
		return false;
//		XYZ n1crd=bs1.getcrd(BackBoneSite::NCRD);
	double theta1 = angle(n1crd, h1crd, o2crd) * 180 / 3.14159265;
	if (!(theta1 > TDHAMIN))return false;
	XYZ c2crd = bs2.getcrd(BackBoneSite::CCRD);
	double anoc = angle(n1crd, o2crd, c2crd) * 180.0 / 3.14159265;
	if (!(anoc > TDABMIN))return false;
	double theta2 = angle(h1crd, o2crd, c2crd) * 180 / 3.14159265;
	if(!(theta2 >THABMIN)) return false;
	if (hbg) {
		hbg->rha = sqrt(dh1o2);
		hbg->rda = sqrt(dn1o2);
		hbg->adha = theta1;
		hbg->ahab = theta2;
		hbg->adab = anoc;
	}
	return true;
}
void NSPproteinrep::segmenthbonded(const std::vector<BackBoneSite> &seg,
		std::vector<unsigned int> &hbonded) {
	hbonded.clear();
	hbonded.resize(seg.size(), 0);
	for (unsigned int i = 0; i < seg.size() - 1; ++i) {
		for (unsigned int j = i + 1; j < seg.size(); ++j) {
			unsigned int hb = NSPproteinrep::hbonded(seg[i], seg[j]);
			if (hb == 1) {
				hbonded[i] += 1;
				hbonded[j] += 10;
			} else if (hb == 2) {
				hbonded[i] += 10;
				hbonded[j] += 1;
			}
		}
	}
}

void NSPproteinrep::updatehbonded(const std::vector<BackBoneSite> &seg1,
		std::vector<unsigned int> &hbonded1,
		const std::vector<BackBoneSite> &seg2,
		std::vector<unsigned int> &hbonded2) {
	for (unsigned int i = 0; i < seg1.size(); ++i) {
		for (unsigned int j = 0; j < seg2.size(); ++j) {
			unsigned int hb = hbonded(seg1[i], seg2[j]);
			if (hb == 1) {
				hbonded1[i] += 1;
				hbonded2[j] += 10;
			} else if (hb == 2) {
				hbonded1[i] += 10;
				hbonded2[j] += 1;
			}
		}
	}
}
