#ifndef NASTRAN_READER_H_
#define NASTRAN_READER_H_

#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cctype>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ostream>
#include "matvec.h"

#define FIELD_LENGTH	8

using namespace std;

void NastranReader ();

// class NastranReader {
// protected:
  // friend class SPH::TriMesh;
  // friend class SPH::Domain;
	// std::vector <std::string> rawData;
	// int line_count;
	// int elem_count;
	// int node_count;
  
  // //Flattened arrays such as GPU type in order of mantain this
  // double  *node;
  // int     *elcon;
	// int 		*nodeid;	//If node number does not begin in one
	// std::map <int,int> nodepos;	//id to position
  
  // //TriMesh trimesh;
	
	// public:
  // int     dim;
  // NastranReader(){dim=3;}
	// NastranReader(char* fName){read(fName);}
	
	// void WriteCSV(char const * FileKey);
	// void WriteVTK(char const * FileKey);
	
  // ~NastranReader();
	// virtual inline void read(char *fName);
	
// };

////void NastranReader::read( char* fName){
  
void ReadNastranTriMesh( char* fName, double *node, int *elcon, int dim, bool issurf);


#endif

