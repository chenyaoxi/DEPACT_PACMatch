/*
 * pdbfile.h
 *
 *  Created on: 2016年2月23日
 *      Author: hyliu
 */

#ifndef PDBFOLDER_H_
#define PDBFOLDER_H_
#include "dataio/datapaths.h"
#include <string>
namespace NSPdataio {
/*!This class handles paths and filenames of pdb files
 * from the protein data bank.The pdb files are organized into subfolders
 * based on (the third letter) of the pdbids.
 *
 */
class PdbFolder {
public:
	static PdbFolder &getinstance(const std::string &root=""){
		static PdbFolder instance(pdbroot(root));
		return instance;
	}
	/*! Checks if a pdb file of the given pdbid exists locally.
	 * It will try to download the from the internet if a localcopy does not exist.
	 * If success, filename with complete path is returned.
	 *
	 */
	std::string getfile(std::string pdbid,
			const std::string & extension="ent");

	/*!Locate path to a pdb file with given subfolder name(the third letter of odbid)
	 * If the subdir does not exists, try to create it, so that later
	 * downloading can be carried out.
	 *
	 */
	std::string subfolder(std::string name);

	/*!root path for pdb files
	 *
	 */
	const std::string & root() const {return root_;}

	std::string & root() {return root_;}
private:
	std::string root_;
	PdbFolder(const std::string & r) {if(r.empty())root_="."; else root_=r;}
};
/*!download a pdb file from the internet with a given local filename (and path)
 * It uses the "wget" program to download the file.
 */
std::string downloadpdb(std::string pdbid, std::string localfile);

}
#endif /* PDBFOLDER_H_ */
