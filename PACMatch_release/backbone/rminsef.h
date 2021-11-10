/*
 * rminsef.h
 *
 *  Created on: 2016年6月8日
 *      Author: hyliu
 */

#ifndef RMINSEF_H_
#define RMINSEF_H_
#include "backbone/backbonesite.h"
#include "backbone/backbonetorsionstate.h"
#include <map>
namespace NSPproteinrep {
class BBD2Matrix {
public:
	enum {
		N = 0, C = 1, CB = 2, CX = 3
	};
	BBD2Matrix(const BackBoneSite &s1, const BackBoneSite &s2);
	BBD2Matrix(const BackBoneSite &s1, const BackBoneSite &s2, BBatomdis2 &atomdis2);
	int min_index() const {
		return min_index_;
	}
	int max_index() const {
		return max_index_;
	}
	std::pair<int, int> minmaxindices() const {
		return std::make_pair(min_index_, max_index_);
	}
	std::vector<double> &matrix() {
		return dist2_;
	}
	double phi() {
		return phi_;
	}
private:
	std::vector<double> dist2_;
	int min_index_;
	int max_index_;
	double phi_;
	std::pair<int, int> indextopair(int idx);
};
class EneTables {
public:
	struct Table {
		std::vector<double> r;
		std::vector<double> ene;
	};
	typedef std::pair<int, int> OrienType;
	typedef std::string SSType;
	typedef std::pair<SSType, OrienType> EneType;
	typedef std::pair<int, int> PairType;
	typedef std::pair<PairType, OrienType> PEType;

	static double energy(const Table &tab, double r);
	static double subenergy(const Table &rtab,
			const std::vector<Table> &ptables, double r, double p);
	/*	double energy(EneType & t, double r) const {
	 if(tables_.find(t) == tables_.end()) return 10.0;
	 return energy(tables_.at(t),r);
	 }
	 */
	double energy(EneType & t, double r, double p) const {
		if (tables_.find(t) == tables_.end())
			return 10.0;
		double etot = energy(tables_.at(t), r);
		if (phiene_)
			etot += subenergy(tables_.at(t), subtables_.at(t), r, p);
		return etot;
	}
	double energy(PEType &t, double r) const {
		if (petables_.find(t) == petables_.end())
			return 10.0;
		return energy(petables_.at(t), r);
	}
	std::map<EneType, Table> & tables() {
		return tables_;
	}
	const std::map<EneType, Table> & tables() const {
		return tables_;
	}
//	void read(const std::string & filename);
	void readwithsubtables(const std::string &filename);
	void readcoiltables(const std::string &filename);
private:
	static int getindex(const std::vector<double> &crds, double c) {
		return intervale(crds, c, 0, crds.size() - 1);
	}
	static int intervale(const std::vector<double> &crds, double c, int lower,
			int upper);
	std::map<EneType, Table> tables_;
	std::map<EneType, std::vector<Table>> subtables_;
	std::map<PEType, Table> petables_;
	bool phiene_ { false };
};

class RMinSEF {
public:
	static RMinSEF &getinstance(const std::string &ssfilename = "",
			const std::string &coilfilename = "") {
		static RMinSEF instance;
		static bool initialized { false };
		if (!initialized) {
			assert(ssfilename != "" && coilfilename != "");
			instance.init(ssfilename);
			instance.initcoiltables(coilfilename);
			initialized=true;
		}
		return instance;
	}
	void init(const std::string & filename) {
		enetables_.readwithsubtables(filename);

	}
	void initcoiltables(const std::string &filename) {
		enetables_.readcoiltables(filename);
	}
	static typename EneTables::EneType enetype(char s1, char s2,
			const BBD2Matrix & d2m);
	static typename EneTables::PEType petype(const BackBoneSite & bs1,
			const BackBoneSite &bs2, const BBD2Matrix & d2m);
	double onebody(BackBoneSite &s) {
		return 0.0;
	}
	double twobody(const BackBoneSite &s1, const BackBoneSite &s2);
//	double twobody(BackBoneSite &s1, BackBoneSite &s2,BBatomdis2 &atomdis2);
	EneTables & enetables() {
		return enetables_;
	}
	const EneTables & enetables() const {
		return enetables_;
	}
private:
	EneTables enetables_;
};
}

#endif /* RMINSEF_H_ */
