/*
 * SasaPSD.cpp
 *
 *  Created on: 2017Äê10ÔÂ23ÈÕ
 *      Author: notxp
 */

#include "designseq/SasaPSD.h"

namespace NSPdesignseq {

SasaPSD::SasaPSD() {
	radii = 3.5;
	psdList[0] = XYZ(-1.040, 0.006, 1.345); /* CB */
	psdList[1] = XYZ(-2.332, -0.663, 0.981); /* CG */
	psdList[2] = XYZ(-3.372, -0.657, 2.325); /* CD */
	points[  0] = XYZ(-1.8054,-1.0317, 5.4323);
	points[  1] = XYZ(-3.1450, 0.3479, 5.6700);
	points[  2] = XYZ(-4.0581,-3.5753, 4.1314);
	points[  3] = XYZ(-2.7891,-1.6531, 5.6293);
	points[  4] = XYZ(-4.8827,-3.3978, 0.7580);
	points[  5] = XYZ(-4.1505, 1.8632, 4.6256);
	points[  6] = XYZ(-5.6259,-1.2039, 4.9463);
	points[  7] = XYZ(-4.9850,-2.9044, 4.4693);
	points[  8] = XYZ(-3.8522,-1.9342, 5.5481);
	points[  9] = XYZ(-2.2757,-3.7019, 3.6580);
	points[ 10] = XYZ(-0.1980,-2.7350, 3.3520);
	points[ 11] = XYZ(-4.6663,-3.8628, 1.7797);
	points[ 12] = XYZ(-4.2060, 0.1091, 5.6368);
	points[ 13] = XYZ(-3.6840,-0.8464, 5.8059);
	points[ 14] = XYZ(-2.5979,-4.0583, 2.6110);
	points[ 15] = XYZ(-3.5992,-4.1473, 2.1965);
	points[ 16] = XYZ(-5.2469, 1.5585, 4.2811);
	points[ 17] = XYZ(-1.8757, 1.7134, 4.4208);
	points[ 18] = XYZ(-5.5415, 2.0582, 1.9114);
	points[ 19] = XYZ(-6.1628, 0.4905, 0.5517);
	points[ 20] = XYZ(-4.1022, 2.5151, 1.0386);
	points[ 21] = XYZ(-3.8115,-3.9185, 1.1337);
	points[ 22] = XYZ(-6.0827,-2.8671, 2.4568);
	points[ 23] = XYZ(-6.6001,-2.0086, 2.2707);
	points[ 24] = XYZ(-5.7765,-2.8931, 3.5367);
	points[ 25] = XYZ(-6.7689,-1.1973, 1.6778);
	points[ 26] = XYZ(-1.8825,-2.1172, 5.1356);
	points[ 27] = XYZ(-5.9851, 1.3715, 1.1817);
	points[ 28] = XYZ(-5.7903,-2.1411, 4.3742);
	points[ 29] = XYZ(-2.0133,-3.0133, 4.5276);
	points[ 30] = XYZ(-6.7875,-1.1418, 2.9163);
	points[ 31] = XYZ(-2.9147,-2.6331, 5.1774);
	points[ 32] = XYZ(-6.5765, 0.5848, 1.6622);
	points[ 33] = XYZ(-2.9764, 2.7795, 2.8578);
	points[ 34] = XYZ(-2.9897, 1.9452, 4.6342);
	points[ 35] = XYZ(-6.6382,-0.2705, 3.5220);
	points[ 36] = XYZ(-6.3656,-1.4768, 0.7075);
	points[ 37] = XYZ(-2.6833,-0.5838, 5.7558);
	points[ 38] = XYZ(-5.7034,-1.5861,-0.1145);
	points[ 39] = XYZ(-4.6766,-1.0809, 5.5450);
	points[ 40] = XYZ(-6.5946, 0.5961, 2.8679);
	points[ 41] = XYZ(-4.4429, 1.9838, 0.0930);
	points[ 42] = XYZ(-2.4358, 1.1429, 5.1770);
	points[ 43] = XYZ(-5.8580, 1.5747, 3.3688);
	points[ 44] = XYZ(-5.1970,-0.1790, 5.2731);
	points[ 45] = XYZ(-3.9835,-2.8616, 4.9738);
	points[ 46] = XYZ(-1.0685,-2.3652, 4.3315);
	points[ 47] = XYZ(-3.9915,-3.6328, 0.1584);
	points[ 48] = XYZ(-4.6504, 2.2667, 3.7630);
	points[ 49] = XYZ(-5.2099, 2.2681, 2.8875);
	points[ 50] = XYZ(-6.0809,-0.2825, 4.5094);
	points[ 51] = XYZ(-1.6128,-3.6551, 2.7334);
	points[ 52] = XYZ(-5.5606, 0.6514, 4.7225);
	points[ 53] = XYZ(-5.3103,-3.5608, 2.5723);
	points[ 54] = XYZ(-3.0415,-3.4320, 4.4323);
	points[ 55] = XYZ(-5.5871,-2.5684, 0.4042);
	points[ 56] = XYZ(-3.5156, 2.8182, 1.9348);
	points[ 57] = XYZ(-1.3272, 0.8072, 4.7591);
	points[ 58] = XYZ(-6.4013,-2.0746, 3.3564);
	points[ 59] = XYZ(-3.5586, 2.4883, 3.8489);
	points[ 60] = XYZ(-3.3322,-3.9820, 3.4173);
	points[ 61] = XYZ(-6.1976, 0.7016, 3.8806);
	points[ 62] = XYZ(-5.6165,-3.2134, 1.5021);
	points[ 63] = XYZ(-1.3005,-3.1419, 3.6607);
	points[ 64] = XYZ(-5.0997, 2.0589, 0.9505);
	points[ 65] = XYZ(-6.8509,-0.2774, 2.2674);
	points[ 66] = XYZ(-5.4044, 1.2660, 0.2223);
	points[ 67] = XYZ(-4.3101,-3.9756, 2.9225);
	points[ 68] = XYZ(-1.1419,-0.2506, 4.9917);
	points[ 69] = XYZ(-5.4815, 0.2733,-0.3083);
	points[ 70] = XYZ(-3.0848,-3.8516,-0.2504);
	points[ 71] = XYZ(-4.6201, 2.5932, 1.9664);
	points[ 72] = XYZ(-4.6406, 0.9833, 5.1446);
	points[ 73] = XYZ(-4.8926,-2.0741, 5.1410);
	points[ 74] = XYZ(-6.0699,-0.5622, 0.0974);
	points[ 75] = XYZ(-6.3670,-1.2330, 4.0421);
	points[ 76] = XYZ(-6.2016, 1.4028, 2.3472);
	points[ 77] = XYZ(-6.2483,-2.3635, 1.2932);
	points[ 78] = XYZ(-1.0361,-1.4168, 4.8182);
	points[ 79] = XYZ(-2.4679, 2.4160, 3.7354);
	points[ 80] = XYZ(-4.1010, 2.7158, 2.9104);
	points[ 81] = XYZ(-2.4256, 3.1770, 1.8696);
	points[ 82] = XYZ(-2.0098, 0.1943, 5.4346);
	points[ 83] = XYZ(-6.6205,-0.3550, 1.0577);
	points[ 84] = XYZ(-4.9536,-3.5105, 3.5923);
	points[ 85] = XYZ(-3.5769, 1.2133, 5.2763);
	points[ 86] = XYZ(-3.4787,-2.0345,-2.0281);
	points[ 87] = XYZ(-3.9532, 1.0095,-1.6314);
	points[ 88] = XYZ(-1.9530,-3.9500,-0.1610);
	points[ 89] = XYZ(-4.2138,-2.4400,-1.3752);
	points[ 90] = XYZ(-2.5450,-3.5305,-1.0145);
	points[ 91] = XYZ(-3.7449,-3.2593,-0.8933);
	points[ 92] = XYZ(-3.6030, 2.2200,-0.5440);
	points[ 93] = XYZ(-4.6056,-0.4822,-1.6738);
	points[ 94] = XYZ(-4.6656,-2.9495,-0.2746);
	points[ 95] = XYZ(-4.7362, 1.2894,-0.6494);
	points[ 96] = XYZ(-3.0558,-2.8742,-1.6338);
	points[ 97] = XYZ(-1.4030,-4.0300, 0.7563);
	points[ 98] = XYZ(-1.7906,-4.0419, 1.7160);
	points[ 99] = XYZ(-0.8490,-3.7160, 1.8370);
	points[100] = XYZ(-4.6680,-1.4490,-1.5050);
	points[101] = XYZ(-2.7210,-4.1250, 0.6410);
	points[102] = XYZ(-5.3370,-0.7750,-0.8100);
	points[103] = XYZ(-4.0160, 0.2710,-1.9410);
	points[104] = XYZ(-3.4180, 2.5770, 0.2250);
	points[105] = XYZ(-5.0340,-2.1170,-0.7030);
	points[106] = XYZ(-2.7600,-4.1020, 1.4750);
	points[107] = XYZ(-4.8850, 0.4100,-1.1590);
	points[108] = XYZ(-3.7840,-1.1030,-2.1730);
	points[109] = XYZ(-3.8440, 1.6550,-1.1610);
	points[110] = XYZ(-0.6550,-3.1770, 2.7460);
	points[111] = XYZ(-1.4360, 2.5550, 3.7110);
	points[112] = XYZ(-1.9340, 3.0620, 2.7980);
	points[113] = XYZ(-0.2000,-0.6110, 4.6860);
	points[114] = XYZ(-0.4880, 1.3660, 4.5220);
	points[115] = XYZ(-0.8810, 2.1060, 4.1410);
	points[116] = XYZ(-0.0340,-1.3690, 4.4020);
	points[117] = XYZ(-0.0570,-2.1060, 3.9570);
	points[118] = XYZ(-0.1560, 0.3330, 4.7160);
	points[119] = XYZ(-2.9780, 2.9000, 0.9940);
	sasaIndex[  0] = 0.101;
	sasaIndex[  1] = 0.228;
	sasaIndex[  2] = 0.272;
	sasaIndex[  3] = 0.303;
	sasaIndex[  4] = 0.327;
	sasaIndex[  5] = 0.347;
	sasaIndex[  6] = 0.365;
	sasaIndex[  7] = 0.379;
	sasaIndex[  8] = 0.393;
	sasaIndex[  9] = 0.405;
	sasaIndex[ 10] = 0.416;
	sasaIndex[ 11] = 0.426;
	sasaIndex[ 12] = 0.435;
	sasaIndex[ 13] = 0.444;
	sasaIndex[ 14] = 0.453;
	sasaIndex[ 15] = 0.461;
	sasaIndex[ 16] = 0.469;
	sasaIndex[ 17] = 0.476;
	sasaIndex[ 18] = 0.483;
	sasaIndex[ 19] = 0.490;
	sasaIndex[ 20] = 0.497;
	sasaIndex[ 21] = 0.504;
	sasaIndex[ 22] = 0.510;
	sasaIndex[ 23] = 0.517;
	sasaIndex[ 24] = 0.523;
	sasaIndex[ 25] = 0.529;
	sasaIndex[ 26] = 0.535;
	sasaIndex[ 27] = 0.541;
	sasaIndex[ 28] = 0.546;
	sasaIndex[ 29] = 0.552;
	sasaIndex[ 30] = 0.558;
	sasaIndex[ 31] = 0.563;
	sasaIndex[ 32] = 0.569;
	sasaIndex[ 33] = 0.574;
	sasaIndex[ 34] = 0.579;
	sasaIndex[ 35] = 0.585;
	sasaIndex[ 36] = 0.590;
	sasaIndex[ 37] = 0.595;
	sasaIndex[ 38] = 0.600;
	sasaIndex[ 39] = 0.605;
	sasaIndex[ 40] = 0.611;
	sasaIndex[ 41] = 0.616;
	sasaIndex[ 42] = 0.621;
	sasaIndex[ 43] = 0.626;
	sasaIndex[ 44] = 0.631;
	sasaIndex[ 45] = 0.636;
	sasaIndex[ 46] = 0.641;
	sasaIndex[ 47] = 0.646;
	sasaIndex[ 48] = 0.651;
	sasaIndex[ 49] = 0.656;
	sasaIndex[ 50] = 0.661;
	sasaIndex[ 51] = 0.666;
	sasaIndex[ 52] = 0.671;
	sasaIndex[ 53] = 0.676;
	sasaIndex[ 54] = 0.680;
	sasaIndex[ 55] = 0.685;
	sasaIndex[ 56] = 0.690;
	sasaIndex[ 57] = 0.695;
	sasaIndex[ 58] = 0.700;
	sasaIndex[ 59] = 0.705;
	sasaIndex[ 60] = 0.709;
	sasaIndex[ 61] = 0.714;
	sasaIndex[ 62] = 0.719;
	sasaIndex[ 63] = 0.724;
	sasaIndex[ 64] = 0.729;
	sasaIndex[ 65] = 0.734;
	sasaIndex[ 66] = 0.739;
	sasaIndex[ 67] = 0.744;
	sasaIndex[ 68] = 0.749;
	sasaIndex[ 69] = 0.754;
	sasaIndex[ 70] = 0.759;
	sasaIndex[ 71] = 0.764;
	sasaIndex[ 72] = 0.769;
	sasaIndex[ 73] = 0.774;
	sasaIndex[ 74] = 0.779;
	sasaIndex[ 75] = 0.784;
	sasaIndex[ 76] = 0.789;
	sasaIndex[ 77] = 0.794;
	sasaIndex[ 78] = 0.799;
	sasaIndex[ 79] = 0.805;
	sasaIndex[ 80] = 0.810;
	sasaIndex[ 81] = 0.815;
	sasaIndex[ 82] = 0.820;
	sasaIndex[ 83] = 0.826;
	sasaIndex[ 84] = 0.831;
	sasaIndex[ 85] = 0.837;
	sasaIndex[ 86] = 0.843;
	sasaIndex[ 87] = 0.849;
	sasaIndex[ 88] = 0.855;
	sasaIndex[ 89] = 0.862;
	sasaIndex[ 90] = 0.869;
	sasaIndex[ 91] = 0.875;
	sasaIndex[ 92] = 0.882;
	sasaIndex[ 93] = 0.889;
	sasaIndex[ 94] = 0.895;
	sasaIndex[ 95] = 0.900;
	sasaIndex[ 96] = 0.905;
	sasaIndex[ 97] = 0.910;
	sasaIndex[ 98] = 0.914;
	sasaIndex[ 99] = 0.919;
	sasaIndex[100] = 0.923;
	sasaIndex[101] = 0.928;
	sasaIndex[102] = 0.932;
	sasaIndex[103] = 0.936;
	sasaIndex[104] = 0.940;
	sasaIndex[105] = 0.944;
	sasaIndex[106] = 0.948;
	sasaIndex[107] = 0.952;
	sasaIndex[108] = 0.955;
	sasaIndex[109] = 0.959;
	sasaIndex[110] = 0.962;
	sasaIndex[111] = 0.966;
	sasaIndex[112] = 0.970;
	sasaIndex[113] = 0.973;
	sasaIndex[114] = 0.977;
	sasaIndex[115] = 0.981;
	sasaIndex[116] = 0.984;
	sasaIndex[117] = 0.988;
	sasaIndex[118] = 0.992;
	sasaIndex[119] = 0.995;
	sasaIndex[120] = 0.998;
}

int SasaPSD::exposeNum(NSPproteinrep::BackBoneSite* resA, vector<NSPproteinrep::BackBoneSite*>& resList){
	vector<XYZ> XYZList;
	LocalFrame cs = getBackboneSiteLocalFrame(*resA);
	XYZ CB = cs.local2globalcrd(psdList[0]);
	XYZ CG = cs.local2globalcrd(psdList[1]);
	XYZ CD = cs.local2globalcrd(psdList[2]);
	float cutoff = 4*radii*radii;
	float radiiSq = radii*radii;
	float dd1, dd2, dd3;
	int resNum = resList.size();
	for(int i=0;i<resNum;i++)
	{
		NSPproteinrep::BackBoneSite* resB = resList.at(i);
		/*
		 * same atom
		 */
		if(resB->cacrd().squaredDistance(cs.origin_) < 0.0001) continue;


		vector<XYZ> resBXYZs;
		resB->getcrd(resBXYZs);

		int bbAtomNum = resBXYZs.size();
		for(int j=0;j<bbAtomNum;j++)
		{
			XYZ& t = resBXYZs.at(j);
			if(CB.squaredDistance(t) > 289) //distance larger than 17 angstrom
				break;
			dd1 = t.squaredDistance(CB);
			dd2 = t.squaredDistance(CG);
			dd3 = t.squaredDistance(CD);
			if(dd1 > cutoff && dd2 > cutoff && dd3 > cutoff)
				continue;
			XYZList.push_back(cs.global2localcrd(t));
		}

		LocalFrame cs2 = getBackboneSiteLocalFrame(*resB);
		for(int j=0;j<3;j++)
		{
			XYZ t = cs2.local2globalcrd(psdList[j]);
			if(CB.squaredDistance(t) > 289) //distance larger than 17 a
				break;
			dd1 = t.squaredDistance(CB);
			dd2 = t.squaredDistance(CG);
			dd3 = t.squaredDistance(CD);
			if(dd1 > cutoff && dd2 > cutoff && dd3 > cutoff)
				continue;
			XYZList.push_back(cs.global2localcrd(t));
		}
	}

	XYZList.push_back(cs.global2localcrd(resA->ocrd()));



	int neighborXYZNum = XYZList.size();

	int n=0;


	for(int i=0;i<120;i++)
	{
		XYZ ball = points[i];
		for(int j=0;j<neighborXYZNum;j++)
		{
			if(XYZList.at(j).squaredDistance(ball) < radiiSq)
			{
				n++;
				break;
			}
		}

	}
	//cout << "point Number: " << neighborXYZNum <<  " expose: " << 120-n <<  endl;
	return 120-n;
}

int SasaPSD::exposeNum(Residue* resA, vector<Residue*>& resList)
{
	vector<XYZ> XYZList;
	LocalFrame cs = resA->getCoordSystem();
	XYZ CB = cs.local2globalcrd(psdList[0]);
	XYZ CG = cs.local2globalcrd(psdList[1]);
	XYZ CD = cs.local2globalcrd(psdList[2]);
	float cutoff = 4*radii*radii;
	float radiiSq = radii*radii;
	float dd1, dd2, dd3;
	int resNum = resList.size();
	for(int i=0;i<resNum;i++)
	{
		Residue* resB = resList.at(i);
		/*
		 * same atom
		 */

		Atom* aca = resA->getAtom("CA");
		if(resB->getAtom("CA")->distance(*aca) < 0.0001) continue;

		vector<XYZ> resBXYZs;
		vector<Atom*>* bbList = resB->getBackboneAtoms();
		for(int j=0;j<bbList->size();j++){
			resBXYZs.push_back(bbList->at(j)->getCoord());
		}

		int bbAtomNum = resBXYZs.size();
		for(int j=0;j<bbAtomNum;j++)
		{
			XYZ& t = resBXYZs.at(j);
			if(CB.squaredDistance(t) > 289) //distance larger than 17 angstrom
				break;
			dd1 = t.squaredDistance(CB);
			dd2 = t.squaredDistance(CG);
			dd3 = t.squaredDistance(CD);
			if(dd1 > cutoff && dd2 > cutoff && dd3 > cutoff)
				continue;
			XYZList.push_back(cs.global2localcrd(t));
		}

		LocalFrame cs2 = resB->getCoordSystem();
		for(int j=0;j<3;j++)
		{
			XYZ t = cs2.local2globalcrd(psdList[j]);
			if(CB.squaredDistance(t) > 289) //distance larger than 17 a
				break;
			dd1 = t.squaredDistance(CB);
			dd2 = t.squaredDistance(CG);
			dd3 = t.squaredDistance(CD);
			if(dd1 > cutoff && dd2 > cutoff && dd3 > cutoff)
				continue;
			XYZList.push_back(cs.global2localcrd(t));
		}
	}

	if(resA->hasAtom("O")){
		XYZList.push_back(cs.global2localcrd(resA->getAtom("O")->getCoord()));
	}


	int neighborXYZNum = XYZList.size();

	int n=0;


	for(int i=0;i<120;i++)
	{
		XYZ ball = points[i];
		for(int j=0;j<neighborXYZNum;j++)
		{
			if(XYZList.at(j).squaredDistance(ball) < radiiSq)
			{
				n++;
				break;
			}
		}

	}
	//cout << "point Number: " << neighborXYZNum <<  " expose: " << 120-n <<  endl;
	return 120-n;
}

int SasaPSD::exposeNumTest(Residue* resA, vector<Residue*>& resList)
{
	cout << "test mode" << endl;
	vector<XYZ> XYZList;
	LocalFrame cs = resA->getCoordSystem();
	XYZ CB = cs.local2globalcrd(psdList[0]);
	XYZ CG = cs.local2globalcrd(psdList[1]);
	XYZ CD = cs.local2globalcrd(psdList[2]);

	printf("TEST: CD %8.3f %8.3f %8.3f\n", CD[0], CD[1], CD[2]);
	printf("TEST: CG %8.3f %8.3f %8.3f\n", CG[0], CG[1], CG[2]);
	printf("TEST: CB %8.3f %8.3f %8.3f\n", CB[0], CB[1], CB[2]);


	float cutoff = 4*radii*radii;
	float radiiSq = radii*radii;
	float dd1, dd2, dd3;
	int resNum = resList.size();
	for(int i=0;i<resNum;i++)
	{
		Residue* resB = resList.at(i);
		/*
		 * same atom
		 */

		Atom* aca = resA->getAtom("CA");
		if(resB->getAtom("CA")->distance(*aca) < 0.0001) continue;

		vector<XYZ> resBXYZs;
		vector<Atom*>* bbList = resB->getBackboneAtoms();
		for(int j=0;j<bbList->size();j++){
			resBXYZs.push_back(bbList->at(j)->getCoord());
		}

		int bbAtomNum = resBXYZs.size();
		for(int j=0;j<bbAtomNum;j++)
		{
			XYZ& t = resBXYZs.at(j);
			if(CB.squaredDistance(t) > 289) //distance larger than 17 angstrom
				break;
			dd1 = t.squaredDistance(CB);
			dd2 = t.squaredDistance(CG);
			dd3 = t.squaredDistance(CD);
			if(dd1 > cutoff && dd2 > cutoff && dd3 > cutoff)
				continue;
			XYZList.push_back(cs.global2localcrd(t));
		}

		LocalFrame cs2 = resB->getCoordSystem();
		for(int j=0;j<3;j++)
		{
			XYZ t = cs2.local2globalcrd(psdList[j]);
			if(CB.squaredDistance(t) > 289) //distance larger than 17 a
				break;
			dd1 = t.squaredDistance(CB);
			dd2 = t.squaredDistance(CG);
			dd3 = t.squaredDistance(CD);
			if(dd1 > cutoff && dd2 > cutoff && dd3 > cutoff)
				continue;
			XYZList.push_back(cs.global2localcrd(t));
		}
	}

	if(resA->hasAtom("O")){
		XYZList.push_back(cs.global2localcrd(resA->getAtom("O")->getCoord()));
	}


	int neighborXYZNum = XYZList.size();

	cout << "point num: " << neighborXYZNum << endl;
	for(int i=0;i<neighborXYZNum;i++){
		printf("point %2d %8.3f %8.3f %8.3f\n", i, XYZList.at(i)[0], XYZList.at(i)[1], XYZList.at(i)[2]);
	}

	int n=0;


	for(int i=0;i<120;i++)
	{
		XYZ ball = points[i];
		for(int j=0;j<neighborXYZNum;j++)
		{
			if(XYZList.at(j).squaredDistance(ball) < radiiSq)
			{
				n++;
				break;
			}
		}

	}
	//cout << "point Number: " << neighborXYZNum <<  " expose: " << 120-n <<  endl;
	return 120-n;
}

SasaPSD::~SasaPSD() {
	// TODO Auto-generated destructor stub
}

} /* namespace NSPdesignseq */
