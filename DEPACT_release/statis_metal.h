/*
 * statis_metal.h
 *
 *  Created on: 2019年12月9日
 *      Author: yxchen
 */

#ifndef STATIS_METAL_H_
#define STATIS_METAL_H_

#include "buildpocket.h"
#include "geometry/xyz.h"
using namespace myobcode;
using namespace subsitedesign;
using namespace OpenBabel;
using namespace NSPproteinrep;
using namespace NSPgeometry;

namespace statismetal
{

bool isresidue(std::string resname);

struct AtmInfo
{
	std::string atmtype;
	XYZ crd;
};
class MetalConf
{
public:
	// read this metal's info
	void readmetal(std::string name, double x, double y, double z)
	{
		metalinfo_.atmtype = name; // residue name
		metalinfo_.crd.x_ = x;
		metalinfo_.crd.y_ = y;
		metalinfo_.crd.z_ = z;
	}
	// add new contactinfo to contactsinfo_
	void addcontact(std::string name, double x, double y, double z)
	{
		AtmInfo ci;
		ci.atmtype = name; // atm name
		ci.crd.x_ = x;
		ci.crd.y_ = y;
		ci.crd.z_ = z;
		contactsinfo_.push_back(ci);
	}
	// output metalinfo
	AtmInfo getmetal() {return metalinfo_;}
	// output the ith contactinfo
	std::vector<AtmInfo> getci() {return contactsinfo_;}
private:
	AtmInfo metalinfo_;
	std::vector<AtmInfo> contactsinfo_;
};
struct ConfStatis
{
	std::vector<double> bonds; // collect all bond lengths: A-Metal; For water, it collect all non-wmbonds.
	std::vector<double> angles; // collect all angles: A-Metal-B
	std::vector<double> wmbonds; // only collect water-metal bond lengths;
//	std::vector<double> diheds; // collect all dihedral angles: Metal-A-B-C
};
std::vector<AAConformersInModel> getconfsfromtitlefile(std::string file, int num_th);
std::vector<AAConformersInModel> getconfsfromtitlefile_pdb(std::string file, int num_th); // read all pdbs
AAConformersInModel getconfsfromtitlefile_spdb(std::string pdbname); // only read single pdb
AAConformersInModel getconfsfromtitlefile_ssdf(std::vector<std::string> sdfname); // only read single sdf
std::vector<MetalConf> aaconf2metalconf(std::vector<AAConformersInModel> cims);
std::vector<MetalConf> aaconf2metalconf_single(AAConformersInModel cim);
std::vector<MetalConf> aaconf2metalconf_spec(std::vector<AAConformersInModel> cims, string metal);
std::vector<MetalConf> aaconf2metalconf_spec_single(AAConformersInModel cim, string metal);
std::map<std::string, ConfStatis> metalconf2statis(std::vector<MetalConf> mcs);
std::map<std::string, ConfStatis> metalconf2statis(std::vector<MetalConf> mcs, std::string spec_metal);
void metalconf2statis_print(std::vector<MetalConf> mcs);
}

#endif /* STATIS_METAL_H_ */
