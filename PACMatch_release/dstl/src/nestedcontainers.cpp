/*
 * nestedcontainers.cpp
 *
 *  Created on: 2015年10月26日
 *      Author: hyliu
 */
#include <cmath>
#include "dstl/nestedcontainers.h"

using namespace NSPdstl;
std::pair<double, double> NSPdstl::columnaverage(
                const NestedVector<double> & table,
                int col) {
    double sum = 0.0;
    double sum2 = 0.0;
    if (table.empty()) return std::make_pair(sum, sum2);
    for (auto iter = table.cbegin(); iter != table.cend();
                    ++iter) {
        double d = (*iter)[col];
        sum += d;
        sum2 += d * d;
    }
    double sz = (double) (table.size());
    return std::make_pair(sum / sz,
                    std::sqrt(
                                    (sum2 - (sum * sum) / sz)
                                                    / sz));
}

