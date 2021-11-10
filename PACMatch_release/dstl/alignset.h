/*
 * alignset.h
 *
 *  Created on: 2017年8月20日
 *      Author: hyliu
 */

#ifndef DSTL_ALIGNSET_H_
#define DSTL_ALIGNSET_H_
#include <vector>
#include <dstl/graph.h>

namespace NSPdstl {

typedef std::vector<std::pair<int,int>> AlignedPositions;
class SetMatch:public Graph{
public:
	virtual bool isneighbor(Vertex v1,Vertex v2) const{
		return !conflict(vertices.at(v1),vertices.at(v2));
	}
	template <typename ELE,typename ABMATCH>
	static AlignedPositions alignset(const std::vector<ELE> &seta,const std::vector<ELE> &setb,
			ABMATCH &abmatch){
		SetMatch sm(seta,setb,abmatch);
		AlignedPositions aln;
		for(auto & ep:maxclique(sm)){
			aln.push_back(sm.indicesab(sm.vertices[ep]));
		}
		return aln;
	};
private:
	template <typename ELE,typename ABMATCH>
	SetMatch(const std::vector<ELE> & seta, const std::vector<ELE> &setb,
			ABMATCH & abmatch){
		sizea_=seta.size();
		sizeb_=setb.size();
		for (int idxa=0;idxa<sizea_;idxa++){
			for (int idxb=0;idxb<sizeb_;idxb++){
				if(abmatch(seta[idxa],setb[idxb])){
					vertices.push_back(makeelementpair(idxa,idxb));
				}
			}
		}
		Graph::size=vertices.size();
	}
	typedef int ElementPair;
	typedef std::vector<ElementPair> Vertices;
	Vertices vertices;
	int sizea_;
	int sizeb_;
	ElementPair makeelementpair(int indexa,int indexb){
		return indexa*sizeb_+indexb;
	}
	std::pair<int,int> indicesab(ElementPair e) const {
		return std::make_pair((int) (e/sizeb_),(int)(e%sizeb_));
	}
	bool conflict(ElementPair ep1,ElementPair ep2) const{
		auto a1b1=indicesab(ep1);
		auto a2b2=indicesab(ep2);
		if(a1b1.first== a2b2.first || a1b1.second == a2b2.second)return true;
		return false;
	}
};


}


#endif /* DSTL_ALIGNSET_H_ */
