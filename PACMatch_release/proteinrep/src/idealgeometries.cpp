/*
 * idealgeometries.cpp
 *
 *  Created on: 2016年11月17日
 *      Author: hyliu
 */
#include "proteinrep/idealgeometries.h"
#include "dataio/inputlines.h"
#include "dataio/splitstring.h"
using namespace NSPproteinrep;
int IdealGeometries::readIdealValues(const std::string &filename) {
	NSPdataio::TextLines lines;
	oldfilename_=filename;
	lines.init(filename);
	int nread = 0;
	for (auto & l : lines.lines())
		if (readLine(l))
			++nread;
	return nread;
}
bool IdealGeometries::readLine(const std::string &line) {
	std::vector<std::string> words = NSPdataio::wordsInString(line);
	if (words.size() == 3)
		return addLength(words);
	if (words.size() == 4)
		return addAngle(words);
	if (words.size() == 5)
		return addRelativeTorsion(words);
	if(words.size()==6)
		return addTorsion(words);
	return false;
}
bool IdealGeometries::addLength(const std::vector<std::string> & words) {
	try {
		double l = std::stod(words[2]);
		Bond b = std::make_pair(words[0], words[1]);
		ideallengths_.insert(std::make_pair(b, l));
		b = std::make_pair(words[1], words[0]);
		ideallengths_.insert(std::make_pair(b, l));
		return true;
	} catch (std::exception &e) {
		return false;
	}
}
bool IdealGeometries::addAngle(const std::vector<std::string> & words) {
	try {
		double val = std::stod(words[3])*degree;
		Angle a(words[0], words[1], words[2]);
		idealangles_.insert(std::make_pair(a, val));
		Angle ar(words[2], words[1], words[0]);
		idealangles_.insert(std::make_pair(ar, val));
		return true;
	} catch (std::exception &e) {
		return false;
	}
}
bool IdealGeometries::addRelativeTorsion(const std::vector<std::string> & words) {
	try {
		double val = std::stod(words[4])*degree;
		RelativeTorsion t(words[0], words[1], words[2], words[3]);
		idealrtorsions_.insert(std::make_pair(t, val));
		RelativeTorsion tr(words[0], words[1], words[3], words[2]);
		idealrtorsions_.insert(std::make_pair(tr, -val));
		return true;
	} catch (std::exception &e) {
		return false;
	}
}
bool IdealGeometries::addTorsion(const std::vector<std::string> & words) {
	try {
		double val = std::stod(words[5])*degree;
		Torsion t(words[1], words[2], words[3], words[4]);
		idealtorsions_.insert(std::make_pair(t, val));
		Torsion tr(words[4], words[3], words[2], words[1]);
		idealtorsions_.insert(std::make_pair(tr, -val));
		return true;
	} catch (std::exception &e) {
		return false;
	}
}

double IdealGeometries::idealLength(const AtomName &ia,
		const AtomName &ja) const {
	Bond b = std::make_pair(ia, ja);
	auto it = ideallengths_.find(b);
	if (it != ideallengths_.end())
		return it->second;
	return -1.0;
}
double IdealGeometries::idealAngle(const AtomName &ia, const AtomName &ja,
		const AtomName &ka) const {
	Angle a(ia, ja, ka);
	auto it = idealangles_.find(a);
	if (it != idealangles_.end())
		return it->second;
	return 1000;
}
double IdealGeometries::idealRelativeTorsion(const AtomName &ja, const AtomName &ka,
		const AtomName &la1, const AtomName &la2) const {
	RelativeTorsion t = RelativeTorsion(ja, ka, la1,la2);
	auto it = idealrtorsions_.find(t);
	if (it != idealrtorsions_.end())
		return it->second;
	return 1000;
}
double IdealGeometries::idealTorsion(const AtomName &ia,const AtomName &ja,
		const AtomName &ka, const AtomName &la) const{
	Torsion t=Torsion(ia,ja,ka,la);
	auto it =idealtorsions_.find(t);
	if(it != idealtorsions_.end()) return it->second;
	return 1000;
}
