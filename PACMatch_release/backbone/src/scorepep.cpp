/*
 * scorepep.cpp
 *
 *  Created on: 2017年6月30日
 *      Author: hyliu
 */

#include "backbone/scorepep.h"
#include "dstl/randomengine.h"
#include "pdbstatistics/phipsidistr.h"
#include "pdbstatistics/proteinblock.h"
using namespace NSPproteinrep;
using namespace domaintree;

static const double D2MAX{625};
static const double D20{25.0};
static const double TREERESOLUTION{5.0};
static const double PGLY{0.1054};
static const double PPRO{0.0623};
static const double NMIN{1.e-8};
static const double RANDOMFRACTION{5.e-3};
double myd2max;
double myd20;
std::vector<std::pair<double,double>> d2maxd0{{
	225,25},{400,64},{900,100},{2025,225},{3600,900}};
static double neighborcut2(int vecsize){
	return myd2max*(double) vecsize;
}
static std::vector<std::string> randomseq(int length){
	std::vector<std::string> seq;
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	rneg.setrealrng(0,1);
	for(int i=0;i<length;++i){
		std::string res="ALA";
		double r=rneg.realrng()();
		if(r<PGLY) res="GLY";
		else if(r<PGLY+PPRO) res="PRO";
		seq.push_back(res);
	}
	return seq;
}
/*
static double refprob_phipsi(double phi,double psi){
	double pgly=NSPpdbstatistics::PhiPsiDistr::glydistr().distr(phi,psi);
	double ppro=NSPpdbstatistics::PhiPsiDistr::transprodistr().distr(phi,psi);
	double pprep=NSPpdbstatistics::PhiPsiDistr::preprodistr().distr(phi,psi);
	double pcoil=NSPpdbstatistics::PhiPsiDistr::coildistr().distr(phi,psi);
	double res=PGLY*pgly+PPRO*pgly+(1-PGLY-PPRO)*
			(PPRO*pprep+(1-PPRO)*pcoil);
	return res;
}
static double refprob_phi(double phi){
	double pgly=NSPpdbstatistics::PhiPsiDistr::glydistr().mdistr_phi(phi);
	double ppro=NSPpdbstatistics::PhiPsiDistr::transprodistr().mdistr_phi(phi);
	double pprep=NSPpdbstatistics::PhiPsiDistr::preprodistr().mdistr_phi(phi);
	double pcoil=NSPpdbstatistics::PhiPsiDistr::coildistr().mdistr_phi(phi);
	double res=PGLY*pgly+PPRO*pgly+(1-PGLY-PPRO)*
			(PPRO*pprep+(1-PPRO)*pcoil);
	return res;
}

static double refprob_psi(double psi){
	double pgly=NSPpdbstatistics::PhiPsiDistr::glydistr().mdistr_psi(psi);
	double ppro=NSPpdbstatistics::PhiPsiDistr::transprodistr().mdistr_psi(psi);
	double pprep=NSPpdbstatistics::PhiPsiDistr::preprodistr().mdistr_psi(psi);
	double pcoil=NSPpdbstatistics::PhiPsiDistr::coildistr().mdistr_psi(psi);
	double res=PGLY*pgly+PPRO*pgly+(1-PGLY-PPRO)*
			(PPRO*pprep+(1-PPRO)*pcoil);
	return res;
}
*/
std::vector<double> ScorePep::extractconf(const std::vector<BackBoneSite> &sites, long posi,int length){
	std::vector<double> result;
	long start=posi-length/2;
	if(start<=0) return result;
	if(!fragstartsite(sites.begin()+start,sites.end(),length,std::string(),false)) return result;
	if(conftype_!=ALL) {
		for(int i=0; i<length-1;++i) {
			double omiga=sites[start+i].omiga();
			if(omiga >-90.0 && omiga <90.0) return result;
		}
	}
	result.push_back(sites[start].psi());
//	result.push_back(sites[start].omiga());
	for(int i=1;i<length-1;++i){
		result.push_back(sites[start+i].phi());
		result.push_back(sites[start+i].psi());
//		result.push_back(sites[start+i].omiga());
	}
	result.push_back(sites[start+length-1].phi());
	for(auto & a:result) {
		while (a <-180.0) a+=360.0;
		while (a>180.0)  a-=360.0;
	}
	char pbtype=NSPpdbstatistics::ProteinBlock::pbtype(result);
	if(conftype_!= ALL){
		if(conftype_==COIL){
			if(pbtype=='m' || pbtype=='d') {
				if (sites[posi].sscodechar() !='C')return std::vector<double>();
			}
		} else if(conftype_ == HELIX &&
			(sites[posi].sscodechar() != 'H' || pbtype !='m')) return std::vector<double>();
		else if (conftype_==STRAND &&
				(sites[posi].sscodechar() !='E' ||pbtype!='d')) return std::vector<double>();
	}
	return result;
}
void ScorePep::buildreftree(long npoints){
	refconfs_=std::shared_ptr<std::vector<std::vector<double>>>
			(new std::vector<std::vector<double>>());
	for (long i=0; i<npoints; ++i) {
		std::vector<double> point=sampleconf(length_);
		refconfs_->push_back(point);
		if(i==0)reftree_.init(refconfs_,TREERESOLUTION);
		reftree_.insertpoint(i);
	}
}
void ScorePep::buildrefphipsitrees(long npoints){
	refweights_=std::shared_ptr<std::vector<double>> (new std::vector<double>());

	refphis_=std::shared_ptr<std::vector<std::vector<double>>>
			(new std::vector<std::vector<double>>());

	refpsis_=std::shared_ptr<std::vector<std::vector<double>>>
				(new std::vector<std::vector<double>>());

	refphipsis_=std::shared_ptr<std::vector<std::vector<double>>>
				(new std::vector<std::vector<double>>());

	refphipsips_=std::shared_ptr<std::vector<std::vector<double>>>
					(new std::vector<std::vector<double>>());

	refphipsins_=std::shared_ptr<std::vector<std::vector<double>>>
					(new std::vector<std::vector<double>>());
	assert(npoints<=templateconfs_->size());
	const typename NSPpdbstatistics::PhiPsiDistr * distr=&(NSPpdbstatistics::PhiPsiDistr::mixcoildistr());
	auto &rng=NSPdstl::RandomEngine<>::getinstance().realrng(0,1);
	for (long i=0; i<npoints; ++i) {
//		std::vector<double> point=templateconfs_->at(i);
		std::vector<double> phi;
		std::vector<double> psi;
		std::vector<double> phipsi(2,0.0);
		std::vector<double> phipsip(2,0.0);
		std::vector<double> phipsin(2,0.0);
//		phi.push_back(point[7]);
		double p,s;
		distr->randomphipsi(rng,&p,&s);
		if(p>180.0)p-=360.0;
		if(s>180.0)s-=360.0;
		phi.push_back(p);
		refphis_->push_back(phi);
		distr->randomphipsi(rng,&p,&s);
		if(p>180.0)p-=360.0;
		if(s>180.0)s-=360.0;
		psi.push_back(s);
		refpsis_->push_back(psi);
		distr->randomphipsi(rng,&p,&s);
		if(p>180.0)p-=360.0;
		if(s>180.0)s-=360.0;
		phipsip[0]=p;
		phipsip[1]=s;
		refphipsips_->push_back(phipsip);
		distr->randomphipsi(rng,&p,&s);
		if(p>180.0)p-=360.0;
		if(s>180.0)s-=360.0;
		phipsi[0]=p;
		phipsi[1]=s;
		refphipsis_->push_back(phipsi);
		distr->randomphipsi(rng,&p,&s);
		if(p>180.0)p-=360.0;
		if(s>180.0)s-=360.0;
		phipsin[0]=p;
		phipsin[1]=s;
		refphipsins_->push_back(phipsin);
		double w=1.0;
//		double w=templateweights_->at(i);
		refweights_->push_back(w);
		refwtot_ +=w;
		if(i==0){
			refphi_.init(refphis_,TREERESOLUTION);
			refpsi_.init(refpsis_,TREERESOLUTION);
			refphipsi_.init(refphipsis_,TREERESOLUTION);
			refphipsip_.init(refphipsips_,TREERESOLUTION);
			refphipsin_.init(refphipsins_,TREERESOLUTION);
		}
		refphi_.insertpoint(i);
		refpsi_.insertpoint(i);
		refphipsi_.insertpoint(i);
		refphipsip_.insertpoint(i);
		refphipsin_.insertpoint(i);
	}
/*	std::ofstream ofs;
	ofs.open("refphi.dat");
	for (auto &p:*refphis_) {
		ofs<< p[0]<< std::endl;
	}
	ofs.close();
	ofs.open("refpsi.dat");
	for (auto &p:*refpsis_) {
		ofs<< p[0]<< std::endl;
	}
	ofs.close();
	ofs.open("refphipsi.dat");
	for (auto &p:*refphipsis_) {
		ofs<< p[0]<<" "<<p[1]<< std::endl;
	}
	ofs.close();
	ofs.open("refphipsip.dat");
	for (auto &p:*refphipsips_) {
		ofs<< p[0]<< " "<<p[1]<<std::endl;
	}
	ofs.close();
	ofs.open("refphipsin.dat");
	for (auto &p:*refphipsins_) {
		ofs<< p[0]<< " "<<p[1]<<std::endl;
	}
	ofs.close();*/
}
void ScorePep::buildtree(const std::vector<BackBoneSite> &sites,int peplength){
	length_=peplength;
	templateconfs_=std::shared_ptr<std::vector<std::vector<double>>>
			(new std::vector<std::vector<double>>());
	templateweights_=std::shared_ptr<std::vector<double>>(new std::vector<double>());
/*	std::vector<double> nhelix(peplength,0.0);
	std::vector<double> ncoil(peplength,0.0);
	std::vector<double> nstrand(peplength,0.0);*/
	for (long posi=peplength/2; posi<sites.size()-peplength/2; ++posi) {
		std::vector<double> point=extractconf(sites,posi,peplength);
		if(point.empty()) continue;
/*		if(!excludedpb_.empty()) {
			char pbtype=NSPpdbstatistics::ProteinBlock::pbtype(point);
			if(excludedpb_.find(pbtype) != excludedpb_.end()) {
				char sscode=sites[posi].sscodechar();
				if(sscode !='C') continue;
			}
		}*/
/*		for(int i=0;i<peplength;++i) {
			char sscode=sites[posi-peplength/2+i].sscodechar();
			if(sscode=='H') nhelix[i] +=1.0;
			else if(sscode=='E') nstrand[i]+=1.0;
			else ncoil[i] +=1.0;
		}*/
		double w=1.0;
		if(conftype_==ALL || conftype_==HELIX){
			char pbtype=NSPpdbstatistics::ProteinBlock::pbtype(point);
			if(pbtype=='m' && pbmweight_<1.0){
				auto &rng=NSPdstl::RandomEngine<>::getinstance();
				rng.setrealrng(0.0,1.0);
				if(rng.randomreal()>pbmweight_) continue;
				w=1.0/pbmweight_;
			}
			if(pbtype=='d' && sites[posi].sscodechar()== 'E') w=0.1;
			if(pbtype=='m' && sites[posi].sscodechar()=='H') w*=0.1;
		}
		templateconfs_->push_back(point);
		templateweights_->push_back(w);
		templatewtot_+=w;
		if(templateconfs_->size()==1)tree_.init(templateconfs_,TREERESOLUTION);
		tree_.insertpoint(templateconfs_->size()-1);
	}
/*	phelix_.resize(peplength,0.0);
	pcoil_.resize(peplength,0.0);
	pstrand_.resize(peplength,0.0);
	for(int i=0;i<5;++i) {
		phelix_[i]=nhelix[i]/(nhelix[i]+ncoil[i]+nstrand[i]);
		pcoil_[i]=ncoil[i]/(nhelix[i]+ncoil[i]+nstrand[i]);
		pstrand_[i]=nstrand[i]/(nhelix[i]+ncoil[i]+nstrand[i]);
	}
	std::cout <<"P_helix:";
	for(auto p:phelix_) std::cout<<"  "<<p;
	std::cout<<std::endl;
	std::cout <<"P_strand:";
	for(auto p:pstrand_) std::cout<<"  "<<p;
	std::cout<<std::endl;
	std::cout <<"P_coil:";
	for(auto p:pcoil_) std::cout<<"  "<<p;
	std::cout<<std::endl;*/
}
double ScorePep::pairscore(const std::vector<double> & conf, const std::vector<double> &tmpl){
	double s=1.0;
	for (int i=0;i<tmpl.size();++i) {
		double si;
		double diff2=conf[i]-tmpl[i];
		while(diff2 <-180.0) diff2+=360.0;
		while(diff2 >180.0) diff2 -=360.0;
		diff2=diff2*diff2;
		if(diff2 <= myd20) si=1.0;
		else if(diff2 >=myd2max) si=0.0;
		else {
//			double x=(diff2-myd20)/(myd2max-myd20);
//			si=1-x*x;
			double x=(diff2-myd2max)/(myd20-myd2max);
			si=x*x;
		}
		s*=si;
	}
	return s;
}
double ScorePep::neighborsum(const std::vector<double> &conf, const std::vector<double> &weights,
		TorsionVectorTree & tree) {
	double bound2=neighborcut2(conf.size());
	std::vector<std::vector<double>> *tvdata=&(tree.gettorsionvectors());
	domaintree::D2Leaf<long,std::vector<std::vector<double>>,AngleCrd>
		d2leaf(tvdata,1000000,bound2);
	tree.gettree().findneighbor(conf,d2leaf,bound2);
	std::vector<std::pair<long,double>> &neighbors=d2leaf.nnearest().neighbors();
	double s=0.0;
//	double s=NMIN;
	for(auto &n:neighbors){
		if(n.second < 0.00001) continue;  //ignore self
		const std::vector<double> & temp=(*tvdata)[n.first];
		s+=pairscore(conf,temp)*weights[n.first];
	}
//	if(s<1.0e-16) s=1.0e-16;
	return s;
}
double ScorePep::refneighborsum(const std::vector<double> &conf) {
	std::vector<double> phi;
	std::vector<double> psi;
	std::vector<double> phipsi(2,0.0);
	std::vector<double> phipsip(2,0.0);
	std::vector<double> phipsin(2,0.0);
	phi.push_back(conf[7]);
	double sumphi=neighborsum(phi,*refweights_,refphi_)/refwtot_;
	psi.push_back(conf[0]);
	double sumpsi=neighborsum(psi,*refweights_,refpsi_)/refwtot_;
	phipsi[0]=conf[3];
	phipsi[1]=conf[4];
	double sumphipsi=neighborsum(phipsi,*refweights_,refphipsi_)/refwtot_;
	phipsip[0]=conf[1];
	phipsip[1]=conf[2];
	double sumphipsip=neighborsum(phipsip,*refweights_,refphipsip_)/refwtot_;
	phipsin[0]=conf[5];
	phipsin[1]=conf[6];
	double sumphipsin=neighborsum(phipsi,*refweights_,refphipsin_)/refwtot_;
	return sumphi*sumpsi*sumphipsi*sumphipsip*sumphipsin;
}
double ScorePep::score(const std::vector<double> &conf,double *sx,double *sr_ana,
		double *d2m,double *d20) {
	double mins=1.e-8;
	double cut=6.0/templatewtot_;
	for(auto& dd:d2maxd0){
		myd2max=dd.first;
		myd20=dd.second;
//	double mins=NMIN/(double) refsize();
//	*srx=neighborsum(conf,reftree_)/(double) refsize();
//	double sr=*srx+mins;
		*sx=neighborsum(conf,*templateweights_,tree_)/templatewtot_;
		*sr_ana=refneighborsum(conf);
	 	 if(*sx>cut || *sr_ana>cut) break;
//	 	*sx=0.0;
//	 	*sr_ana=0.0;
//	double npseudo=1.0/templatesize();
//	double s=(*sx+(*sr_ana)*npseudo)/((double)(templatesize())+npseudo);
//	*sx/=(double) (templatesize());
	}
	double s=*sx+mins;
	double sr=*sr_ana+mins;
	*d2m=myd2max;
	*d20=myd20;
	return s/sr;
}
std::vector<double> ScorePep::sampleconf(int length,int peptype){
	std::vector<double> point;
	std::vector<std::string> rseq=randomseq(length+1);
	NSPdstl::RandomEngine<> &rneg = NSPdstl::RandomEngine<>::getinstance();
	rneg.setrealrng(0,1);
/*	long imax=refphis_->size();
		long ir=rneg.intrng(0,imax)();
		point.push_back(refpsis_->at(ir)[0]);
		ir=rneg.intrng()();
		point.push_back(refphipsips_->at(ir)[0]);
		point.push_back(refphipsips_->at(ir)[1]);
		ir=rneg.intrng()();
		point.push_back(refphipsis_->at(ir)[0]);
		point.push_back(refphipsis_->at(ir)[1]);
		ir=rneg.intrng()();
		point.push_back(refphipsins_->at(ir)[0]);
		point.push_back(refphipsins_->at(ir)[1]);
		ir=rneg.intrng()();
		point.push_back(refphis_->at(ir)[0]);*/
	 double p=rneg.randomreal();
	 const NSPpdbstatistics::PhiPsiDistr *distr;
	 if((p<1.0/3.0 || peptype==HELIX) && peptype!=STRAND
			 &&peptype!=COIL) distr=&(NSPpdbstatistics::PhiPsiDistr::helixdistr());
	 else if((p<2.0/3.0 || peptype== STRAND) && peptype!=COIL) distr=&(NSPpdbstatistics::PhiPsiDistr::stranddistr());
	 else distr=&(NSPpdbstatistics::PhiPsiDistr::mixcoildistr());
	 for(int i=0;i<length;++i){
		double phi,psi;
		distr->randomphipsi(rneg.realrng(), &phi, &psi);
		if(i==0) point.push_back(psi);
		else if(i==length-1) point.push_back(phi);
		else {
			point.push_back(phi);
			point.push_back(psi);
		}
	 }
	for(auto &a:point) {
		while(a<-180.0) a+=360.0;
		while(a>180.0) a-=360.0;
	}
	return point;
}
/*
double ScorePep::refprobability(const std::vector<double> &conf) const {
	double sum=0;
	std::vector<const NSPpdbstatistics::PhiPsiDistr *> distrs;
	distrs.push_back(&(NSPpdbstatistics::PhiPsiDistr::glydistr()));
	distrs.push_back(&(NSPpdbstatistics::PhiPsiDistr::transprodistr()));
	distrs.push_back(&(NSPpdbstatistics::PhiPsiDistr::preprodistr()));
	distrs.push_back(&(NSPpdbstatistics::PhiPsiDistr::coildistr()));
	std::vector<double> pseq;
	pseq.push_back(PGLY);
	pseq.push_back(PPRO);
	pseq.push_back((1-PGLY-PPRO)*PPRO);
	pseq.push_back(1-(pseq[0]+pseq[1]+pseq[2]));
	double ptotal=1.0;
	for(int i=0;i<length_;++i){
		double p=0.0;
		if(i==0) {
			for(int d=0;d<4;++d)
				p += pseq[d]*distrs[d]->mdistr_psi(conf[0]);
		} else if(i==length_-1) {
			for(int d=0; d<4;++d)
				p +=pseq[d]*distrs[d]->mdistr_phi(conf.back());
		} else {
			for(int d=0;d<4; ++d)
				p +=pseq[d]*distrs[d]->distr(conf[2*i-1],conf[2*i]);
		}
//		std::cout <<p<<" ";
		ptotal *=p;
	}
//	std::cout <<std::endl;
	return ptotal;
}*/
