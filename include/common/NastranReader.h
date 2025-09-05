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
//#include "matvec.h"

#define FIELD_LENGTH	8
#include <iomanip>

using namespace std;

namespace MetFEM {


using namespace std;

class Domain_d;  
class TriMesh_d;
class NastranReader {
protected:
  friend class TriMesh_d;
  friend class Domain_d;
	std::vector <std::string> rawData;
	int line_count;
	int elem_count;
	int node_count;
  
  //Flattened arrays such as GPU type in order of mantain this
  double  *node = nullptr;
  int     *elcon= nullptr;
	int 		*nodeid= nullptr;	//If node number does not begin in one
	std::map <int,int> nodepos;	//id to position
  
  //TriMesh trimesh;
	
	public:
  int     dim;
  NastranReader(){dim=3;}
	NastranReader(const char* fName){read(fName);}
	
	void WriteCSV(char const * FileKey);
	void WriteVTK(char const * FileKey);
	
  ~NastranReader();
	virtual inline void read(const char *fName);
	
};

};

#endif

