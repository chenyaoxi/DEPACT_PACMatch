/*
 * symmatrix1d.h
 *
 *  Created on: 2017年8月17日
 *      Author: hyliu
 */

#ifndef DSTL_SYMMATRIX1D_H_
#define DSTL_SYMMATRIX1D_H_

namespace NSPdstl{
template<typename DATA>
class SymMatrix1D: public std::vector<DATA>  {
public:
	static int idx1d (int i,int j) {
		int col=i<j?i:j;
		int row=i<j?j:i;
		return (row+1)*row/2+col;
	}
	static std::pair<int,int> rowcol(int idx){
		int row=0;
		while(idx >=(row+1)*row/2) {
			++row;
		}
		row -=1;
		int col=idx-(row+1)*row/2;
		return std::make_pair(row,col);
	}
	DATA & operator()(int i,int j) {
		return this->at(idx1d(i,j));
	}
	const DATA & operator()(int i,int j) const  {
		return this->at(idx1d(i,j));
	}
	DATA sum() const {
		DATA sum{0.0};
		for(auto &d:*this){
			sum+=d;
		}
		return sum;
	}
	void init(int size){this->clear();this->resize((size+1)*size/2);}
};
}


#endif /* DSTL_SYMMATRIX1D_H_ */
