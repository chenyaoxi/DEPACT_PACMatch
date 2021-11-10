/*
 * readpdb.h
 *
 *  Created on: 2016年11月16日
 *      Author: hyliu
 */

#ifndef PROTEINREP_PDBREADER_H_
#define PROTEINREP_PDBREADER_H_
#include "dstl/mapkeyint.h"
#include "proteinrep/pdbrecord.h"
#include <cassert>
#include <map>
#include <memory>
namespace NSPproteinrep {
class MapPdbKeyInt;
class PdbReader {
public:
	typedef std::pair<int,char> ResKeyType; //residue number + insertion code
	typedef std::map<ResKeyType,std::vector<PdbRecord>> ResMapType;
	typedef std::map<char,ResMapType> RecordsType;
		// atom records organized by chain and residue
	void readpdb(const std::string & filename);
	void readpdb(std::vector<std::string> & lines);
    RecordsType & records() {return records_;}
    const std::vector<std::string> &getaminoacidsequence(char chainid) const{
    	return aminoacidsequences_.at(chainid);
    }
    std::string chainids() const {
    	std::string res;
    	for(auto & c:records_) res.push_back(c.first);
    	return res;
    }
    const RecordsType & records() const {return records_;}
    void addRecord(const PdbRecord & record);
    std::shared_ptr<MapPdbKeyInt> mappdbkeyint() const {return mappdbkeyint_;}
    bool isfailure() {return failure_;}
private:
	RecordsType records_;
	std::map<char,std::vector<std::string>> aminoacidsequences_;
	std::shared_ptr<MapPdbKeyInt> mappdbkeyint_;
	bool failure_ {false};
};

class MapPdbKeyInt {
public:
	MapPdbKeyInt(const typename PdbReader::RecordsType & records);
	char pdbChainID(unsigned int chainnumber=0)
		{return mapchainidint_.key(chainnumber);}
	typename PdbReader::ResKeyType pdbResKey(unsigned int posinumber,unsigned int chainnumber=0) {
		return mapreskeyint_.at(chainnumber).key(posinumber);
	}
	int pdbResID(unsigned int posinumber, unsigned int chainnumber=0){
		return pdbResKey(posinumber, chainnumber).first;
	}
	unsigned int chainNumber(char pdbchainid){
		return mapchainidint_.keynumber(pdbchainid);
	}
	unsigned int posiNumber(const PdbReader::ResKeyType &reskey, char pdbchainid=' ') {
		if (chainNumber(pdbchainid) == -10)
			return -10;
		return mapreskeyint_.at(chainNumber(pdbchainid)).keynumber(reskey);
	}
	const NSPdstl::MapKeyInt<char> & mapchainidint() const {return mapchainidint_;}
	const NSPdstl::MapKeyInt<typename PdbReader::ResKeyType> mapreskeyint(unsigned int chainNumber) const {
		return mapreskeyint_.at(chainNumber);
	}
private:
	NSPdstl::MapKeyInt<char> mapchainidint_;
	std::vector<NSPdstl::MapKeyInt<typename PdbReader::ResKeyType>> mapreskeyint_;
};

struct PDBModelID {
	std::string pdbid{""};
	int biounit{1};
	int model{0};
};


std::string extractpdbmodel(const PDBModelID & mid);
std::string extractpdbmodel(const std::string & filename,const PDBModelID &mid);
std::shared_ptr<PdbReader> readpdbmodel(const std::string & pdbid,int biounit=1,
		int modelno=0);
}

#endif /* PROTEINREP_PDBREADER_H_ */
