/*
 * backbonegeometry.h
 *
 *  Created on: 2016年11月17日
 *      Author: hyliu
 */

#ifndef PROTEINREP_IDEALGEOMETRIES_H_
#define PROTEINREP_IDEALGEOMETRIES_H_
#include "dataio/datapaths.h"
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <cassert>
namespace NSPproteinrep {

class IdealGeometries {
public:
	const double degree{3.14159265358979323846/180.0};
	IdealGeometries() {;}
	static IdealGeometries & getGlobalInstance(const std::string &filename = "idealgeometries.dat") {
		static IdealGeometries *instance_{nullptr};
		if( !instance_) {
#ifdef _OPENMP
#pragma omp critical(idealgeometry)
			{
				if(!instance_){
#endif
			instance_=new IdealGeometries;
//			std::string datapath = NSPdataio::datapath();
//			std::string completefilename=datapath+filename;
			std::string completefilename=NSPdataio::datafilename(filename);
			instance_->readIdealValues(completefilename);
#ifdef _OPENMP
#pragma omp flush
				}
			}
#endif
		}
		if( filename != "idealgeometries.dat" && filename != instance_->oldfilename_) {
			instance_->readIdealValues(filename);
		}
		return *instance_;
	}
	typedef std::string AtomName;
	typedef std::pair<AtomName, AtomName> Bond;
	struct Angle {
		Angle(const AtomName & i, const AtomName & j, const AtomName & k) :
				ia(i), ja(j), ka(k) {
			;
		}
		AtomName ia;
		AtomName ja;
		AtomName ka;
		friend bool operator <(const Angle &a1, const Angle &a2){
			return std::tie(a1.ia,a1.ja,a1.ka) < std::tie(a2.ia,a2.ja,a2.ka);
		}
	};
	struct RelativeTorsion {
		RelativeTorsion(const AtomName & j, const AtomName &k, const AtomName &l1,
				const AtomName &l2) :
				ja(j), ka(k), la1(l1), la2(l2) {
			;
		}
		AtomName ja;
		AtomName ka;
		AtomName la1;
		AtomName la2;
		friend bool operator <(const RelativeTorsion &r1, const RelativeTorsion &r2){
				return std::tie(r1.ja,r1.ka,r1.la1,r1.la2)
					<std::tie(r2.ja,r2.ka,r2.la1,r2.la2);
			}
	};
	struct Torsion {
		Torsion(const AtomName & i, const AtomName &j, const AtomName &k,
				const AtomName &l) :
				ia(i), ja(j), ka(k), la(l) {;}
		AtomName ia;
		AtomName ja;
		AtomName ka;
		AtomName la;
		friend bool operator <(const Torsion &r1, const Torsion &r2){
				return std::tie(r1.ia,r1.ja,r1.ka,r1.la)
					<std::tie(r2.ia,r2.ja,r2.ka,r2.la);
			}
	};
	int readIdealValues(const std::string &filename);
	bool readLine(const std::string & line);
	bool addLength(const std::vector<std::string> &words);
	bool addAngle(const std::vector<std::string> &words);
	bool addRelativeTorsion(const std::vector<std::string> &words);
	bool addTorsion(const std::vector<std::string>  &words);
	double idealLength(const AtomName &ia, const AtomName &ja) const;
	double idealAngle(const AtomName &ia, const AtomName &ja,
			const AtomName &ka) const;
	double idealTorsion(const AtomName &ia,const AtomName &ja,
			const AtomName &ka, const AtomName &la) const;
	double idealRelativeTorsion(const AtomName &ja, const AtomName &ka,
			const AtomName &la1, const AtomName &la2) const;
	std::string filename() const {return oldfilename_;}
private:
	IdealGeometries(const IdealGeometries & g){;}
	std::string oldfilename_ { "" };
	std::map<Bond, double> ideallengths_;
	std::map<Angle, double> idealangles_;
	std::map<RelativeTorsion, double> idealrtorsions_;
	std::map<Torsion,double> idealtorsions_;
};
}

#endif /* PROTEINREP_IDEALGEOMETRIES_H_ */
