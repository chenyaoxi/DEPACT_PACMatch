/*
 * splitstring.h
 *
 *  Created on: 2016年11月4日
 *      Author: hyliu
 */

#ifndef DATAIO_SPLITSTRING_H_
#define DATAIO_SPLITSTRING_H_
#include <string>
#include <vector>
#include <ctype.h>
namespace NSPdataio{

/*!This function locates boundaries (or starting positions of
 * substrings) in a string using a user-supplied cutter.
 * It returns a list of starting positions.
 * if cutter(str[i],str[i+1]) is true, i+1 would be included in the
 * returned list as the starting position of a substring.
 *
 */
template <typename ISCUT>
std::vector<int> splitString(const std::string & str, const ISCUT & cutter){
	std::vector<int> substarts;
	substarts.push_back(0);
	for(int i=0; i<str.length()-1;++i){
		if(cutter(str[i],str[i+1])) substarts.push_back(i+1);
	}
	return substarts;
}

/*!cutter functor
 * Cut at interfaces between digital and non-digital characters
 */
class DigitCutter{
public:
	bool operator()(char a,char b) const {
		return !(isdigit(a) == isdigit(b));
	}
};
/*!cutter functor
 * Cut at alphabet and non-alphabet interfaces
 */
class AlphaCutter{
public:
	bool operator()(char a,char b) const {
		return !(isalpha(a) == isalpha(b));
	}
};
/*!cutter functor
 * Cut at delimiter-non-delimiter interfaces
 */
class DelimCutter{
public:
	DelimCutter(const std::string & del=""):delim_(del){;}
	bool isDelim(char a) const {return isspace(a) || delim_.find(a) != std::string::npos;}
	bool operator()(char a,char b) const {
		return !(isDelim(a) == isDelim(b));
	}
private:
	std::string delim_;
};
/*!extract integers contained in a string
 *
 */
std::vector<int> integersInString (const std::string &str);

/*!
 * extract alphabetic words in a string
 */
std::vector<std::string> alphaWordsInString (const std::string &str);

/*!
 * extracts words separated by one or more delimiting characters in a string
 */
std::vector<std::string> wordsInString (const std::string &str,const std::string &delim="");
}

#endif /* DATAIO_SPLITSTRING_H_ */
