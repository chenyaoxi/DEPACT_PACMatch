/*
 * testinputlines.cpp
 *
 *  Created on: 2015年10月21日
 *      Author: hyliu
 */
#include <iostream>
#include "dstl/nestedcontainers.h"
#include "dataio/inputlines.h"

int main(int argc, char **argv) {
    using namespace NSPdataio;
    std::string str="5, 8, 9, -12, 16-33";
      std::vector <int> selected=stringtoselection(str);
      for(auto i:selected) std::cout <<"\t"<<i;
      std::cout<<std::endl;
 /*   InputLines il;
    std::string inputfile("../../testdata/testin.txt");
    if (argc >= 2) inputfile = std::string(argv[1]);
//    il.init(inputfile, '#', std::vector<int>{3,5,5,5});
    il.init(inputfile, '#');
    std::unique_ptr<InputLines> psublines = il.extractLines(
                    "EP");*/
    /*	for (auto line=psublines->begin(); line !=psublines->end(); line++) {
     for (auto wrd=line->begin(); wrd != line->end(); wrd++) {
     std::cout << *wrd <<" ";
     }
     std::cout <<std::endl;
     }
     */
/*    il.printElements();
    il.removeLines("EP");
    il.printElements();
    NSPdstl::NestedVector<double> table =
                    il.extractTable();
    table.printElements();*/
}

