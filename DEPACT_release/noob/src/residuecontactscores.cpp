/*
 * tmpltresiduescores.cpp
 *
 *  Created on: 2018年9月4日
 *      Author: hyliu
 */
#include "residuecontactscores.h"
#include "atomtypes.h"
#include "scorecontact.h"

using namespace subsitedesign;

ResidueContactScores::ResidueContactScores(const TmpltContacts & tcons,
		const std::vector<std::string> &atypes){
	const std::vector<AtomContacts> &acons=*(tcons.atomcontacts);
	const std::vector<std::vector<NSPproteinrep::AAConformer>> &aas =
					tcons.model->conformers;
	setup(acons,aas,atypes);
}
void ResidueContactScores::setup(const std::vector<AtomContacts> &acons,
		const std::vector<std::vector<NSPproteinrep::AAConformer>> &aas,
		const std::vector<std::string> &atypes)
{
	int nlatom=atypes.size();
	for(int i=0; i<acons.size();++i){
		const AtomContacts &ac=acons.at(i);
		for(auto &d:ac.dcontacts){
			ResidueID rid(d.ichain,d.iresidue);
			const NSPproteinrep::AAConformer &aa=
					aas.at(d.ichain).at(d.iresidue);
			std::string rname=aa.residuename;
			for(auto &ad2:d.iatomd2){
				std::string paname=rname;
				std::string patype=rname;
				if(aa.isnaturalaa()){
					paname=aa.atomlist.at(ad2.first.second);
					patype=proteinatomtype(rname,paname);
				}
				std::string latype=atypes.at(ad2.first.first);
				double r=sqrt(ad2.second);
				double s=ScoreContact::score(latype,patype,r);
				if(aa.isnaturalaa())
					insertplscore(rid,paname,nlatom,ad2.first.first,s);
				else insertmlscore(rid,nlatom,ad2.first.first,s);
			}
		} //dcontacts
		for(auto &md:ac.mcontacts){
			ResidueID mid=md.first;
			std::string matype=aas.at(mid.first).at(mid.second).residuename;
			for(auto &d:md.second){
				ResidueID rid(d.ichain,d.iresidue);
				const NSPproteinrep::AAConformer &aa=
						aas.at(d.ichain).at(d.iresidue);
				std::string rname=aa.residuename;
				for(auto &ad2:d.iatomd2){
					std::string paname=aa.
								atomlist.at(ad2.first.second);
					std::string patype=proteinatomtype(rname,paname);
					double r=sqrt(ad2.second);
					double score=ScoreContact::score(matype,patype,r);
					insertmpscore(mid,rid,paname,score);
				}
			}
		}//mcontacts
	} //acons
}
