/*
 * StringTool.h
 *
 *  Created on: 2017��10��17��
 *      Author: notxp
 */

#ifndef DESIGNSEQ_STRINGTOOL_H_
#define DESIGNSEQ_STRINGTOOL_H_

#include <string>
#include <vector>

namespace NSPdesignseq {

using namespace std;

string trimString(const string& s);
void splitString(const string& s, string delim, vector<string >* ret);

}
#endif /* DESIGNSEQ_STRINGTOOL_H_ */
