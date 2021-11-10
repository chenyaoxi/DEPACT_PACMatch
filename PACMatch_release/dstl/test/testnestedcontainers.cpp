/*
 * testnestedcontainer.cpp
 *
 *  Created on: 2016年11月4日
 *      Author: hyliu
 */

#include "dstl/nestedcontainers.h"

int main() {
    using namespace NSPdstl;
    NestedVector<double> nvd;
    std::vector<double> l1{5,4,3,2};
    std::vector<double> l2{0,1.2};
    nvd.push_back(std::move(l1));
    nvd.push_back(std::move(l2));
    nvd.printElements();
}


