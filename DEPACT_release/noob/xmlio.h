/*
 * xmlio.h
 *
 *  Created on: 2018年8月9日
 *      Author: hyliu
 */

#ifndef XMLIO_H_
#define XMLIO_H_
#include "tmpltssas.h"
#include <libxml++/libxml++.h>
#include <sstream>
#include <cassert>
namespace subsitedesign {
/**
 * @collection of fucntions storing and retrieving substructure alignments as XML files
 *
 * @TODO probably better to use boost archiving functions
 */
TmpltSSAs::Alignment xml2alignment(xmlpp::Node *ele) ;
Move3D xml2move3d(xmlpp::Node *ele);
TmpltSSAs xml2tmpltssas(xmlpp::Node *e);
TargetStruct xml2targetstruct(xmlpp::Node *e);
void addnodemove3d(xmlpp::Node *parent,const Move3D &move3d);
void addnodealignment(xmlpp::Node *parent, const TmpltSSAs::Alignment &algn);
void addnodetmpltssas(xmlpp::Node *parent, const TmpltSSAs &ssas);
}

#endif /* XMLIO_H_ */
