/*
 * inputlines.cpp
 *
 *  Created on: 2015年9月28日
 *      Author: hyliu
 */
/*! Implements InputLines
 *
 */
#include "dataio/inputlines.h"
#include <iostream>
#include <fstream>
#include <utility>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm>

#define MAXLINELENGTH 1024  //maximum length of one line in the input file

using namespace NSPdataio;
std::pair<std::string,std::string> NSPdataio::keyvalstrings(const std::string &line){
  auto res=std::make_pair<std::string,std::string>("","");
    int key_len=line.find_first_of('=');
    if(key_len == std::string::npos) return res;
    res.first=boost::trim_copy(line.substr(0,key_len));
    res.second=boost::trim_copy(line.substr(key_len+1));
    return res;
}
typedef boost::tokenizer<boost::char_separator<char> > Tokenizer_char;
typedef boost::tokenizer<boost::offset_separator> Tokenizer_int;
std::vector<std::string> NSPdataio::parseline(
                const std::string & inputline,
                const std::vector<int> & cols) {
    std::vector<std::string> parsedline;
    if (cols.empty()) {
        boost::char_separator<char> sep(" ,\t;");/*! char delimiter set is  " ,\t;"*/
        Tokenizer_char tok(inputline, sep);
        for (auto it = tok.begin(); it != tok.end(); ++it)
            parsedline.push_back(std::string(*it));

    } else {
        boost::offset_separator offsets(cols.begin(),
                        cols.end());
        Tokenizer_int tok(inputline, offsets);
        for (auto it = tok.begin(); it != tok.end(); ++it)
            parsedline.push_back(
                            boost::trim_copy(
                                            std::string(
                                                            *it)));
    }
    return parsedline;
}

std::vector<int> NSPdataio::stringtoselection(const std::string &str){
	std::vector<int> res;
	std::vector<std::string> words=parseline(str,std::vector<int>());
	try {
		for(auto & w:words){
			int posi=w.find('-');
			if(posi==std::string::npos) {
				res.push_back(std::stoi(w));
			} else {
				int rend=std::stoi(w.substr(posi+1));
				int rstart;
				if(posi==0) {
					rstart=res.back()+1;
				} else {
					rstart=std::stoi(w.substr(0,posi));
				}
				for(int i=rstart;i<=rend;++i){
					res.push_back(i);
				}
			}
		}
	}catch (std::exception &e) {
			std::cout<<"Error processing selection string in stringtoselection.The input string is" <<std::endl;
			std::cout <<"\""<<str<<"\""<<std::endl;
			return std::vector<int>();
	}
	std::unique(res.begin(),res.end());
	return res;
}

void NSPdataio::attendError(const std::string &errormessage) {
    std::cout << errormessage << std::endl;
    exit(1);        //exit on error here
}

void TextLines::init(const std::string &filename) {
    std::ifstream ifs;
    resize(0);
    char inputline[MAXLINELENGTH];
    ifs.open(filename.c_str(), std::ifstream::in);
    if (ifs.fail()) {
        std::string msg = "Failed to open inputfile "
                        + filename;
        attendError(msg);
    }
    while (ifs.getline(inputline, MAXLINELENGTH)) {
        push_back(std::string(inputline));
        if (back().empty()) pop_back();        //do not store empty lines
    }
    ifs.close();
}
;

void InputLines::init(const std::string &filename,
                char commentchar,
                const std::vector<int> & cols) {
    std::ifstream ifs;
    char inputline[MAXLINELENGTH];
    resize(0);
    ifs.open(filename.c_str(), std::ifstream::in);
    if (ifs.fail()) {
        std::string msg = "Failed to open inputfile "
                        + filename;
        attendError(msg);
    }
    while (ifs.getline(inputline, MAXLINELENGTH)) {
        if (inputline[0] == commentchar) continue;
        push_back(parseline(std::string(inputline), cols));
        if (back().empty()) pop_back();        //do not store empty lines
    }

    if (not ifs.eof()) {
        std::string msg = "Error occurred reading "
                        + filename;
        attendError(msg);
    }
    ifs.close();
}

std::unique_ptr<InputLines> InputLines::extractLines(
                const std::string & label) const {
    std::unique_ptr<InputLines> psublines(new InputLines);
    for (auto lineit = cbegin(); lineit != cend();
                    lineit++) {
        bool keep = false;
        for (auto sit = lineit->begin();
                        sit != lineit->end(); sit++) {
            if ((*sit) == label) {
                keep = true;
                break;
            }
        }
        if (keep) psublines->push_back(*lineit);
    }
    return psublines;
}

void InputLines::removeLines(const std::string & label) {
    for (auto lineit = begin(); lineit != end();) {
        bool remove = false;
        for (auto sit = lineit->begin();
                        sit != lineit->end(); sit++) {
            if (*sit == label) {
                remove = true;
                break;
            }
        }
        if (remove) lineit = nestedvector().erase(lineit);
        else lineit++;
    }
}

NSPdstl::NestedVector<double> InputLines::extractTable() const {
    NSPdstl::NestedVector<double> table;
    typedef convertelements<std::string, double> CE;
    CE ce(&table);
    auto fun = static_cast<std::function<
                    void(const std::string &)> >(std::bind(
                    &CE::convertOneEle, ref(ce),
                    std::placeholders::_1));
    auto nl = static_cast<std::function<void()> >(std::bind(
                    &CE::lineEnd, ref(ce)));
    processElements(fun, nl);
    /*
     for (auto it=begin(); it!=end();t++){
     std::vector<double> array= wordsToNumbers<double>(*it);
     if(array.size() > 0) {
     table.push_back(std::move(array));
     }
     }
     */
    return table;
}

