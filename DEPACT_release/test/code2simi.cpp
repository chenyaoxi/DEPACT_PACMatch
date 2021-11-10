/*
 * code2simi.cpp
 *
 *  Created on: 2020年7月30日
 *      Author: yxchen
 */

#include "noob/contactcode.h"
#include "buildpocket.h"
#include "noob/theozyme.h"
//extern "C" {
//#include "cliquer/cliquer.h"
//#include "cliquer/graph.h"
//}

using namespace NSPproteinrep;
using namespace ContactCode;
using namespace std;

/*
 * calc. similarity for many codes.
 * input codeA.txt, codeB.txt angle(degree)
 * output << codeA.txt codeB.txt similarity <<
 */

int main(int argc, char **argv)
{
// read all info.
	if (argc <=3)
	{
		std::cout << "input codeA.txt codeB.txt distance_threshold(Angstrom)" << std::endl;
		exit(1);
	}
	double d_th = stod(argv[3]);
	string fileA = argv[1];
	string fileB = argv[2];
	PocCode pcA, pcB;
	pcA.readcode(fileA);
	pcB.readcode(fileB);

// calc. simi by Tanimoto.
	// calc. intersection.
	auto cA = pcA.contacts();
	auto cB = pcB.contacts();
	auto lcA = pcA.ligcrds();
	auto lcB = pcB.ligcrds();
	auto ecA = pcA.envcrds();
	auto ecB = pcB.envcrds();
	std::set<int> matchA, matchB;
	for (int pa = 0; pa < cA.size(); pa++)
		for (int pb = 0; pb < cB.size(); pb++)
		{
			if (lcA[cA[pa].i].x_ != lcB[cB[pb].i].x_
					|| lcA[cA[pa].i].y_ != lcB[cB[pb].i].y_
					|| lcA[cA[pa].i].z_ != lcB[cB[pb].i].z_)
				continue;
			if (ecA[cA[pa].j].distance(ecB[cB[pb].j]) > d_th) continue;
			bool simi = false;
			for (auto ta : cA[pa].types)
				for (auto tb : cB[pb].types)
					if (ta == tb) simi = true;
					else if (ta == 3 && tb < 3) simi = true;
					else if (tb == 3 && ta < 3) simi = true;
			if (simi)
			{
				matchA.insert(pa);
				matchB.insert(pb);
			}
		} // compare every pa & pb
	int insertion = min(matchA.size(), matchB.size()); // = min(cA_match, cB_match)
	// calc. union.
	int uni = cA.size() + cB.size() - insertion;
	// similarity & coverage
	std::cout << fileA << " " << fileB << " " << (double)insertion/uni << " " << (double)insertion/cA.size() << std::endl;
}


