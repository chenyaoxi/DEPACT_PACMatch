/*
 * datapaths.h
 *
 *  Created on: 2017年4月25日
 *      Author: hyliu
 */

#ifndef DATAIO_DATAPATHS_H_
#define DATAIO_DATAPATHS_H_
#include <string>
namespace NSPdataio{
/*!get path string from a system ENVIRONMENTAL variable
 * returns an empty string if envvar is not a defined environmental variable name.
 */
std::string getenvpath(const std::string & envvar);

/*!get datapath using the environmental variable "ABACUS_DATAPATH"
 *
 */
std::string datapath();

/*! if filename exists in the current working directory, returns filename
 *  otherwise return datapath()+filename
 */
std::string datafilename(const std::string &filename);

/*!get path storing downloaded pdb files using the environmental variable
 * "DOWNLOADED_PDB_PATH"
 *
 */
std::string downloadedpdbpath();
/*!
 * get root path of pdb files
 */
std::string pdbroot(const std::string &root);
}


#endif /* DATAIO_DATAPATHS_H_ */
