/*
 * backbonealignment.h
 *
 *  Created on: 2017年8月20日
 *      Author: hyliu
 */
#include "backbone/backbonesite.h"
#include "geometry/rigidtransform.h"
#include "dstl/alignset.h"

#include <set>
namespace NSPproteinrep {

inline std::pair<int,int> dualindex(int singleindex,int size){
	return std::make_pair<int,int>(singleindex/size,singleindex%size);
}
inline int singleindex(int index1,int index2,int size){
	return index1*size+index2;
}

class BackBoneAlignment{
public:
	enum {AUTO_SEEDS,SIMPLE_SEEDS,SS_SEEDS};
	BackBoneAlignment(int seedmode=AUTO_SEEDS):seedmode_(seedmode){;}
	struct Results{
		std::shared_ptr<NSPgeometry::RigidTransform> rigidtransform;
		NSPdstl::AlignedPositions alignedpositions;
		double rmsd2;
		int length() const {return alignedpositions.size();}
		static bool better(const Results & res_ref, const Results &res){
			if(res.rmsd2<=9.0) {
				if(res.length()==res_ref.length()){
					return res.rmsd2<res_ref.rmsd2;
				} else
					return res.length()>res_ref.length();
			}
			return false;
		}
	};
	void init(const std::vector<BackBoneSite> &confa, const std::vector<BackBoneSite> &confb);
	static std::shared_ptr<Results> align(const std::vector<BackBoneSite> &confa,
			const std::vector<BackBoneSite> &confb,int seedmode={AUTO_SEEDS});
	void optalign(Results *operes);
	void selectseeds();
	void simple_selectseeds();
	void SS_selectseeds(const std::vector<int> &ssseqa,const std::vector<int> &ssseqb);
	void tryalign(int seeda,int seedsb,Results *res);
private:
	const std::vector<BackBoneSite> *confa_{nullptr};
	const std::vector<BackBoneSite> *confb_{nullptr};
	std::vector<NSPgeometry::XYZ> crda_;
	std::vector<NSPgeometry::XYZ> crdb_;
	std::vector<int> seedsa_;
	std::vector<int> seedsb_;
	std::vector<std::set<int>> treatedseedb_;
	void updatetreated(const std::vector<std::pair<int,int>> &ap);
	int seedmode_;
};
NSPdstl::AlignedPositions alignsites(const std::vector<BackBoneSite> &seta, const std::vector<BackBoneSite> &setb,
		double rcut);
}



