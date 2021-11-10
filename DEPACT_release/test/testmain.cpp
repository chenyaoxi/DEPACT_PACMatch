/*
 * testmain.cpp
 *
 *  Created on: 2018年8月6日
 *      Author: hyliu
 */


#include <iostream>

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
using namespace OpenBabel;
int main(int argc,char **argv)
{
  OBMol mol;
  OBConversion obconversion;
  obconversion.SetInFormat("sdf");
  std::string file("ipt.sdf");
  bool notatend = obconversion.ReadFile(&mol,file);
  while (notatend)
  {
    std::cout << "Molecular Weight: " << mol.GetMolWt() << std::endl;

    mol.Clear();
    notatend = obconversion.Read(&mol);
  }

  return(0);
}

