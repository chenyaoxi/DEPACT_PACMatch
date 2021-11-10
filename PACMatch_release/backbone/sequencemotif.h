/*
 * sequencemotif.h
 *
 *  Created on: 2016年12月28日
 *      Author: hyliu
 */

#ifndef BACKBONE_SEQUENCEMOTIF_H_
#define BACKBONE_SEQUENCEMOTIF_H_

namespace NSPproteinrep{

class SequenceMotif {
public:
	template<typename ITER>
	static std::string tomotif(ITER iter,int length){
		std::string seqmotif;
		std::string omigaseq;
		for(int i=0;i<length; ++i){
			if(iter->resname=="GLY") motif=motif+"G";
			else if(iter->resname=="PRO"){
				motif=motif+"P";
			}
			else motif=motif+"X";
			if(iter->omiga() <-90|| iter->omiga()>90){
				omigaseq +="c";
			} else {
				omigaseq +="t";
			}
			++iter;
		}
		return seqmotif+omigaseq;
	}
};

}



#endif /* BACKBONE_SEQUENCEMOTIF_H_ */
