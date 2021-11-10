/*
 * minimizer.h
 *
 *  Created on: 2017年8月11日
 *      Author: hyliu
 */

#ifndef DSTL_MINIMIZER_H_
#define DSTL_MINIMIZER_H_
#include <memory>
#include <cassert>
#include <vector>
namespace NSPdstl {
struct MinimizerControl{
	enum{CONSTANTTEMP, TRIPHASE};
	int maxsteps{-1};
	double temperature{-1.0};
	int maxnochangesteps{-1};
	int nprint{1000};
	double temp_high{10.0};
	double temp_low{0.1};
	std::vector<int> phaselengths;
	void settriphase(double th,double tl,const std::vector<int> & pl){
		assert(pl.size()==3);
		temp_high=th;
		temp_low=tl;
		phaselengths=pl;
		mode=TRIPHASE;
	}
	int mode{CONSTANTTEMP};
	double gettemperature(int step) const {
		if(mode==CONSTANTTEMP)  return temperature;
		else if(mode ==TRIPHASE){
			step = step%phaselengths[2];
			if(step<phaselengths[0]) return temp_high;
			else if(step<phaselengths[1]) {
				double alpha = (double)(step-phaselengths[0])/(double)(phaselengths[1]-phaselengths[0]);
				return alpha *temp_low +(1.-alpha)*temp_high;
			}else
				return temp_low;
		}
	}
	bool stop(int nsteps,int nochangesteps) const {
		return nsteps >=maxsteps || nochangesteps >= maxnochangesteps;
	};
};
template<typename OBJ>
class SA_Minimizer{
public:
	virtual std::shared_ptr<OBJ> run(const MinimizerControl &control,const OBJ & obj)=0;
	virtual ~SA_Minimizer(){;}
	double emin() const {return emin_;}
	double e0() const {return e0_;}
protected:
	int nstep_{0};
	double temperature_;
	double emin_{1.e20};
	double e0_{1.e20};
	int nochangesteps_{0};
	std::shared_ptr<OBJ> minconf_;
};

}


#endif /* DSTL_MINIMIZER_H_ */
