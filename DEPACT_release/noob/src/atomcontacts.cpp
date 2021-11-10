/*
 * subsite.cpp
 *
 *  Created on: 2018年8月10日
 *      Author: hyliu
 */

#include "atomcontacts.h"
#include "atomtypes.h"
#include "tmpltssas.h"
#include "scorecontact.h"
#include <fstream>
using namespace subsitedesign;
using namespace NSPproteinrep;
#define DIST2CUTPROTEIN 25.0
#define DIST2CUTNONPROTEIN 25.0

void TmpltContacts::readpdb(bool calccontacts){
	std::vector<std::string> parsedname= TmpltSSAs::parsetmpltname(
			std::string(tmpltname));
	int modelno=std::atoi(parsedname[TmpltSSAs::MODELNO].c_str());
	char chainid=parsedname[TmpltSSAs::CHAINID][0];
	int residuenumber=std::atoi(parsedname[TmpltSSAs::RESIDUENO].c_str());
	auto reader=NSPproteinrep::readpdbmodel(parsedname[TmpltSSAs::PDBID],1,modelno);
	if(!reader) {
		std::cout <<"readpdb failed for " << tmpltname <<std::endl;
		model=nullptr;
		atomcontacts=nullptr;
		return;
	}
	NSPproteinrep::MapPdbKeyInt & mki=*(reader->mappdbkeyint());
	try{
		int ichain=mki.chainNumber(chainid);
		//reskey is pair of resiude sequence and insertion code from PDB record
		int iresidue=mki.posiNumber(std::pair<int,char>(residuenumber,' '),chainid);
		model=std::shared_ptr<NSPproteinrep::AAConformersInModel> (new AAConformersInModel());
		model->getconformers(*reader);
		if(calccontacts)atomcontacts=findcontacts(
			model.get(),ichain,iresidue);
	} catch (std::exception &e){
		model=nullptr;
		atomcontacts=nullptr;
	}
}
std::shared_ptr<std::vector<AtomContacts>> subsitedesign::findcontacts(
		const AAConformersInModel *mdl,int ichain,int iresidue){
	std::shared_ptr<std::vector<AtomContacts>> res=
			std::shared_ptr<std::vector<AtomContacts>>(
					new std::vector<AtomContacts>());
	const AAConformer  & focusconf=mdl->conformers.at(ichain).at(iresidue);
	for (int a=0;a<focusconf.atomlist.size();++a){
		res->push_back(AtomContacts(mdl,ichain,iresidue,a));
	}
	return res;
}

std::shared_ptr<std::vector<AtomContacts>> subsitedesign::findcontacts(
		const AAConformersInModel *mdl,int ichain,int iresidue,
		std::vector<int> selectedatms){
	std::shared_ptr<std::vector<AtomContacts>> res=
			std::shared_ptr<std::vector<AtomContacts>>(
					new std::vector<AtomContacts>());
	const AAConformer  & focusconf=mdl->conformers.at(ichain).at(iresidue);
//	for (int a=0;a<focusconf.atomlist.size();++a){
	for (int a = 0; a < selectedatms.size(); ++a)
	{
		res->push_back(AtomContacts(mdl,ichain,iresidue,a));
	}
	return res;
}


AtomContacts subsitedesign::mergeatomcontacts(const std::vector<AtomContacts> & acontacts,
		const std::vector<int> & atoms){
	AtomContacts res;
	res.model=acontacts.at(0).model;
	for(int a:atoms){
		if (a >= acontacts.size()) {
			std::cout << "al exceeds acontacts' size" << std::endl;
			continue;
		}
		auto & dcontacts=acontacts.at(a).dcontacts;
		for(auto &c:dcontacts){
			bool included=false;
			for( auto &ce:res.dcontacts){
				if(c.sameresidue(ce)){
					for(auto & atomc:c.iatomd2){
						ce.iatomd2.insert(atomc);
					}
					included=true;
					break;
				}
			} //existing contacts
			if(!included){
				res.dcontacts.push_back(c);
			}
		} // dcontacts;
		for(auto &c:acontacts.at(a).mcontacts){
			if(res.mcontacts.find(c.first) == res.mcontacts.end())
				res.mcontacts.insert(c);
		}
	} //atoms
	return res;
}
AtomContacts::AtomContacts(const AAConformersInModel *mdl,
		int focuschain,int focusresidue,
		int focusatom):model(mdl){
	const std::vector<std::vector<AAConformer>> & chains=
			mdl->conformers;
	const AAConformer  & focusconf=chains.at(focuschain).at(focusresidue);
	NSPgeometry::XYZ xyz0=focusconf.globalcrd.at(focusconf.atomlist.at(focusatom));
	for(int c=0;c<chains.size();++c){
		const std::vector<AAConformer> & chain=chains.at(c);
		for(int r=0;r<chain.size();++r){
			if(c==focuschain && r==focusresidue) continue;
			const AAConformer &aa=chain.at(r);
			if(!aa.isnaturalaa() &&aa.atomlist.size()>1) continue; //ignore multi-atom-non protein moities
			double dcut2=aa.isnaturalaa()?DIST2CUTPROTEIN:DIST2CUTNONPROTEIN;
			bool contact=false;
			std::map<std::pair<int,int>,double> contactatoms;
			for(int i=0;i<aa.atomlist.size();++i){
				if(aa.atomlist.at(i)[0]=='H' && aa.isnaturalaa()) continue; // ignore protein H atoms
				NSPgeometry::XYZ xyz1=aa.globalcrd.at(aa.atomlist.at(i));
				double dist2=(xyz0-xyz1).squarednorm();
				if(dist2>dcut2) continue;
				contact=true;
				contactatoms[std::make_pair(focusatom,i)]=dist2;
			}
			if(contact){
				Contact nc(focuschain,focusresidue); //direct contact
				nc.ichain=c;
				nc.iresidue=r;
				nc.iatomd2=contactatoms;
				dcontacts.push_back(nc);
				if(!aa.isnaturalaa()){ //add indirect contact
					NSPgeometry::XYZ xyzw=aa.globalcrd.at(aa.atomlist[0]);
					std::vector<Contact> mcontact;
					for(int mc=0;mc<chains.size();++mc){
						const std::vector<AAConformer> & mchain=chains.at(mc);
						for(int mr=0;mr<mchain.size();++mr){
							const AAConformer &maa=mchain.at(mr);

							// test m-m: water-metal complex should be kept
							if(!maa.isnaturalaa())
							{
								if (maa.residuename != "HOH") continue;
								NSPgeometry::XYZ xyz1=maa.globalcrd.at(maa.atomlist.at(0));
								double dist2=(xyzw-xyz1).squarednorm();
								if(dist2>6.25 || dist2==0) continue; // only accept HOH-metal < 2.5 angstrom
							}

							std::map<std::pair<int,int>,double> mcontactatoms;
							bool mcon=false;
							for(int i=0;i<maa.atomlist.size();++i){
								if(maa.atomlist.at(i)[0]=='H') continue; //ignor protein H atoms
								NSPgeometry::XYZ xyz1=maa.globalcrd.at(maa.atomlist.at(i));
								double dist2=(xyzw-xyz1).squarednorm();
								if(dist2>DIST2CUTNONPROTEIN) continue;
								mcon=true;
								mcontactatoms[std::make_pair(0,i)]=dist2;
							}
							if(mcon){
								Contact ncon(c,r);
								ncon.ichain=mc;
								ncon.iresidue=mr;
								ncon.iatomd2=mcontactatoms;
								mcontact.push_back(ncon);
							}
						} //mr
					}//mc
					if(!mcontact.empty()){
						mcontacts[std::make_pair(c,r)]=mcontact;
					}
					else
						dcontacts.erase(dcontacts.end()); // only media, no p, should be erase
				}//end nonnatural;
			}//end if contact
		} //r
	}//c
}
std::set<int> subsitedesign::bondedligandatoms(const std::vector<ContactDetails> &dtls){
	std::set<int> result;
	for (auto &d:dtls){
		if(d.indirect!=0) continue;
		if(d.distance<1.9) result.insert(d.laid);
	}
	return result;
}


std::set<int> subsitedesign::extbondedligandatoms(const std::vector<ContactDetails> &dtls,
		const std::vector<std::vector<int>> &con){
			std::set<int> batoms=bondedligandatoms(dtls);
// first and second bonded neighbors of these atoms will also be ignored
			if(!batoms.empty()){
				for(int d=0;d<2;++d){
					std::set<int> bck=batoms;
					for(auto &s:bck){
						for(auto c:con[s]){
							batoms.insert(c);
						}
					}
				}
			}
			return batoms;
		}
bool keepcontact(const ContactDetails &d) {
	static std::set<std::string> mediators { "HOH","MN","MG","ZN","CA" };
//	static std::set<std::string> mediators { "HOH"};
	static std::set<std::string> ignoretypes{"E_other","X_eother","Un0"};
	if(ignoretypes.find(d.latype) != ignoretypes.end()) return false;
	if(d.indirect!=0){
		if (mediators.find(d.mligand) == mediators.end())	return false;
	}
	return true;
}
std::vector<ContactDetails> subsitedesign::collectdetails(const TmpltContacts &tc,
		const std::vector<std::vector<std::string>> &atypes){
	std::vector<std::string> parsedname = TmpltSSAs::parsetmpltname(
			tc.tmpltname);
	ContactDetails details;
	std::vector<ContactDetails> dtls;
	details.pdbid = parsedname[TmpltSSAs::PDBID];
	for (int i = 0; i < atypes.size(); ++i) {
		if (atypes[i].size() != 1)
			continue;
		std::string atype = atypes.at(i).at(0);
		if (AtomType::ishydrogen(atype))
			continue;
		NSPproteinrep::AAConformersInModel & mdl = *(tc.model);
		std::vector<std::vector<NSPproteinrep::AAConformer>> &aas =
				mdl.conformers;
		AtomContacts & cts = tc.atomcontacts->at(i);
		std::map<std::pair<int, int>, double> mdistances;
		for (auto &c : cts.dcontacts) {
			details.lcid = c.lchain;
			details.lrid = c.lresidue;
			details.laid = i;
			details.latype = atype;
			details.indirect = 0;
			NSPproteinrep::AAConformer &aa = aas.at(c.ichain).at(
					c.iresidue);
			if (!(aa.isnaturalaa())) {
				mdistances[std::make_pair(c.ichain, c.iresidue)] = sqrt(
						c.iatomd2.begin()->second);
				continue;
			}
			details.pcid = c.ichain;
			details.prid = c.iresidue;
			details.prname = aa.residuename;
			for (auto& aac : c.iatomd2) {
				details.paname = aa.atomlist[aac.first.second];
				details.distance = sqrt(aac.second);
//					ofs<<details.tostring()<<std::endl;
//					ofs.flush();
				dtls.push_back(details);
			}
		}
		for (auto &mc : cts.mcontacts) {
			details.indirect = 1;
			details.mcid = mc.first.first;
			details.mrid = mc.first.second;
			details.mligand =
					aas.at(mc.first.first).at(mc.first.second).residuename;
			details.mdistance = mdistances.at(mc.first);
			for (auto &c : mc.second) {
				NSPproteinrep::AAConformer &aa = aas.at(c.ichain).at(
						c.iresidue);
				details.pcid = c.ichain;
				details.prid = c.iresidue;
				details.prname = aa.residuename;
				for (auto &aac : c.iatomd2) {
					details.paname = aa.atomlist[aac.first.second];
					details.distance = sqrt(aac.second);
//						ofs<<details.tostring()<<std::endl;
//						ofs.flush();
					dtls.push_back(details);
				}
			}
		} //indrect
	} //atypes
	return dtls;
}
double subsitedesign::scorelpcontact(const ContactDetails &dtl){
	std::string latype=dtl.latype;
	if(dtl.indirect!=0) latype=dtl.mligand;
	std::string patype=proteinatomtype(dtl.prname,dtl.paname);
	double r=dtl.distance;
	return ScoreContact::score(latype,patype,r);
}
double subsitedesign::scorelmcontact(const ContactDetails &dtl){
	assert(dtl.indirect!=0);
	std::string latype=dtl.latype;
	std::string matype=dtl.mligand;
	double r=dtl.mdistance;
	return ScoreContact::score(latype,matype,r);
}
std::vector<ContactDetails> subsitedesign::getdtlsnextmol(std::ifstream &ifs,
		std::string &molname,int &nlatoms,
		std::vector<std::string> &atypes,std::set<int> &excludedatoms) {
	static int linenumber{0};
	std::vector<ContactDetails> dts;
	nlatoms = 0;
	bool nextmol = true;
	bool success=true;
	while (ifs.good()) {
		std::string line;
		std::getline(ifs, line);
		++linenumber;
		if (!ifs)
			break;
		if (nextmol) {
			std::stringstream lsstr(line);
			int atsize;
			lsstr >> molname >> nlatoms >> atsize;
			atypes.resize(atsize);
			for (int i = 0; i < atsize; ++i)
				lsstr >> atypes[i];
			nextmol = false;
			std::string nline;
			std::getline(ifs,nline);
			++linenumber;
			lsstr.clear();
			lsstr.str(nline);
			int nexcl;
			lsstr >>nexcl;
			for(int i=0;i<nexcl;++i){
				int at;
				lsstr>>at;
				excludedatoms.insert(at);
			}
			dts.clear();
			continue;
		}
		if (line.substr(0, 6) == "ENDMOL") {
			break;
		} else {
			ContactDetails details;
			std::istringstream iss(line);
			details.read(iss);
/*			if(iss.fail()) {
				std::cout<<"Error reading contact details at line " <<linenumber<<std::endl;
				exit(1);
			}*/
			if(iss.fail()) success=false;
			else if (keepcontact(details))
				dts.push_back(details);
		}
	}
	if(!success) dts.clear();
	return dts;
}
