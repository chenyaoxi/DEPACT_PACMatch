/*
 * randomengine.h
 *
 *  Created on: 2016年4月21日
 *      Author: hyliu
 */

#ifndef RANDOMENGINE_H_
#define RANDOMENGINE_H_
#include <memory>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace NSPdstl {
#ifdef _OPENMP
template <class Engine,class Distribution>
class Generator_Safe: public boost::variate_generator<Engine,Distribution>{
public:
	typedef boost::variate_generator<Engine,Distribution> Generator;
	Generator_Safe(Engine e, Distribution d,omp_lock_t *l):
		Generator(e,d),lock_(l){;}
	typename Distribution::result_type operator()(){
		typename Distribution::result_type result;
		omp_set_lock(lock_);
		result=Generator::operator()();
		omp_unset_lock(lock_);
		return result;
	}
private:
	omp_lock_t *lock_;
};
#endif

template<typename RNGENGINE = boost::mt19937>
struct RandomEngine {
public:
#ifdef _OPENMP
	typedef Generator_Safe<RNGENGINE&, boost::uniform_real<double>> real_generator_type;
	typedef Generator_Safe<RNGENGINE&, boost::uniform_int<long>> int_generator_type;
	typedef Generator_Safe<RNGENGINE&, boost::normal_distribution<double>> normal_generator_type;
#else
	typedef boost::variate_generator<RNGENGINE&, boost::uniform_real<double>> real_generator_type;
	typedef boost::variate_generator<RNGENGINE&, boost::uniform_int<long>> int_generator_type;
	typedef boost::variate_generator<RNGENGINE&, boost::normal_distribution<double>> normal_generator_type;
#endif
	static RandomEngine &getinstance() {
		static RandomEngine instance;
		if (!instance.initialized_) {
#ifdef _OPENMP
#pragma omp critical(randomengine_global)
			{
				if(!instance.initialized_){
					instance.init();
#pragma omp flush
				}
			}
#else
			instance.init();
#endif
		}
		return instance;
	}
	template<typename T>
	void init(const T& seed) {
		init();
		reseed(seed);
	}
	void init() {
		engine_ = std::shared_ptr < RNGENGINE > (new RNGENGINE());
		setrealrng(0.0, 1.0);
#ifdef _OPENMP
		omp_init_lock(&lock_);
#endif
		initialized_=true;
	}
	template<typename T>
	void reseed(const T & sd) {
#ifdef _OPENMP
		omp_set_lock(&lock_);
#endif
		engine_->seed(sd);
#ifdef _OPENMP
		omp_unset_lock(&lock_);
#endif
	}
	template<typename T1, typename T2>
	void setrealrng(T1 rmin, T2 rmax) {
		boost::uniform_real<double> pdf1(rmin, rmax);
		realrng_.reset();
#ifdef _OPENMP
		realrng_ = std::shared_ptr < real_generator_type
				> (new real_generator_type(*engine_, pdf1,&lock_));
#else
		realrng_ = std::shared_ptr < real_generator_type
				> (new real_generator_type(*engine_, pdf1));
#endif
	}

	template<typename T1, typename T2>
	void setintrng(T1 imin, T2 imax) {
		boost::uniform_int<long> pdf(imin, imax);
		intrng_.reset();
#ifdef _OPENMP
		intrng_ = std::shared_ptr < int_generator_type
				> (new int_generator_type(*engine_, pdf,&lock_));
#else
		intrng_ = std::shared_ptr < int_generator_type
				> (new int_generator_type(*engine_, pdf));
#endif
	}
	double randomreal() {
		return (*realrng_)();
	}
	long randomint() {
		return (*intrng_)();
	}
	int_generator_type & intrng() {return *intrng_;}
	template<typename T1, typename T2>
	int_generator_type & intrng(T1 imin,T2 imax) {
		setintrng(imin,imax);
		return *intrng_;}


	real_generator_type & realrng() {return *realrng_;}
	template<typename T1, typename T2>
	real_generator_type & realrng(T1 rmin,T2 rmax) {
		setrealrng(rmin,rmax);
		return *realrng_;}
	normal_generator_type & normalrng() {
		if(!normalrng_) {
			boost::normal_distribution<double> nd;
#ifdef _OPENMP
		normalrng_=std::shared_ptr<normal_generator_type>
						(new normal_generator_type(*engine_,nd,&lock_));
#else
			normalrng_=std::shared_ptr<normal_generator_type>
			(new normal_generator_type(*engine_,nd));
#endif
		}
		return *normalrng_;
	}
	double randomnormal(double mean=0.0,double sd=1.0){
		return mean+sd*normalrng()();
	}
	std::shared_ptr<RNGENGINE> engine_;
	std::shared_ptr<real_generator_type> realrng_;
	std::shared_ptr<int_generator_type> intrng_;
	std::shared_ptr<normal_generator_type> normalrng_{nullptr};
#ifdef _OPENMP
	omp_lock_t lock_;
#endif
	RandomEngine() {
		;
	}
	RandomEngine(RandomEngine const&) = delete;
	void operator=(RandomEngine const&) = delete;
private:
	bool initialized_{false};
};
}
#endif /* RANDOMENGINE_H_ */
