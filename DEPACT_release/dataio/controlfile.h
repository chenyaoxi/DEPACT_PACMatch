/*
 * controlfile.h
 *
 *  Created on: 2017年8月19日
 *      Author: hyliu
 */

#ifndef DATAIO_CONTROLFILE_H_
#define DATAIO_CONTROLFILE_H_
#include "dstl/globalinstances.h"
#include <vector>
#include <map>
#include <set>
#include <bitset>
#include <string>
#include <iostream>
#include <fstream>
namespace NSPdataio{
class ControlFile;
typedef NSPdstl::GlobalInstances<ControlFile> ControlFileInstances;

/*!The class stores lines collected from a file.
 * The file usually contains lines specifying "key=value" controlling parameters.
 * The lines are partitioned into blocks.
 * Each block has a name identifier, and stores the lines in one std::vector object.
 *
 * By using the "START blockname" and "END blockname" directives,
 * one control file may contain several blocks.
 * It may also include other control files through the "INCLUDE filename" directive.
 *
 * This class can be used to process controlling input.
 * Then the blocks of control lines can be used to construct objects of
 * Controls or TypedControls objects.
 */
class ControlFile {
public:
	/*!status flags
	 */
	enum {readundefined,queryundefined,filefailure};

	/*!Returns the control lines with a given block name.
	 *
	 */
	std::vector<std::string> getcontrolines(const std::string & controlname) const {
		if(definedcontrolnames_.find(controlname)== definedcontrolnames_.end())
			status_.set(queryundefined);
		if(controllines_.find(controlname) != controllines_.end())
			return controllines_.at(controlname);
		else
			return std::vector<std::string>();
	}
	/*!
	 * check status
	 */
	bool readundefinedname() const {return status_[readundefined];}
	bool queryundefinedname()const {return status_[queryundefined];}
	bool openfilefailed() const {return status_[filefailure];}

	/*!Find out block names (read from an input file) that are not predefined (in the program)
	 */
	std::vector<std::string> undefinednamesread() const {
		std::vector<std::string> results;
		for(auto & m:controllines_){
			if(definedcontrolnames_.find(m.first) == definedcontrolnames_.end())
				results.push_back(m.first);
		}
		return results;
	}
	/*!read control lines from an input stream. The control lines are organized into blocks
	 * in the input file.
	 *
	 */
	void readfile(std::istream &is);

	/*!opens a control file and read its content.
	 * If anything wrong, the "filefailure" bit of the status_ flag will be set.
	 */
	void readfile(const std::string &filename){
		std::ifstream ifs;
		ifs.open(filename.c_str());
		if(!ifs.good()) {
			status_.set(filefailure);
			return;
		}
		readfile(ifs);
		ifs.close();
	}
	/*!add a controlname (block name) to the predefined names.
	 *
	 */
	void addcontrolname(const std::string &controlname){
		definedcontrolnames_.insert(controlname);
	}

	/*!write the control lines to a output stream
	 *
	 */
	void write(std::ostream &os) const;

private:
	/*!Control lines stored in the object
	 *
	 */
	std::map<std::string,std::vector<std::string>> controllines_;
	/*!predefined control names (block names)
	 *
	 */
	std::set<std::string> definedcontrolnames_;
	mutable std::bitset<8> status_;
};

}



#endif /* DATAIO_CONTROLFILE_H_ */
