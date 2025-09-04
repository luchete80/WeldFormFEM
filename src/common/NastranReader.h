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


namespace SPH {


using namespace std;

class Domain;  
class TriMesh;
class NastranReader {
protected:
  friend class SPH::TriMesh;
  friend class SPH::Domain;
	std::vector <std::string> rawData;
	int line_count;
	int elem_count;
	int node_count;
  
  //Flattened arrays such as GPU type in order of mantain this
  double  *node;
  int     *elcon;
	int 		*nodeid;	//If node number does not begin in one
	std::map <int,int> nodepos;	//id to position
  
  //TriMesh trimesh;
	
	public:
  int     dim;
  NastranReader(){dim=3;}
	NastranReader(char* fName){read(fName);}
	
	void WriteCSV(char const * FileKey);
	void WriteVTK(char const * FileKey);
	
  ~NastranReader();
	virtual inline void read(char *fName);
	
};

void NastranReader::read( char* fName){
	string fileName = fName;
  string line;
  //rawData="";
	fstream file;
    bool found=false;
	//MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	cout << "[I] Reading ... "<<endl;
	file.open(fileName.c_str());
	if (file.is_open()) {
		//cout << "[I] Found input file " << fileName << endl;
		found=true;
	} else {
		//cerr << "[E] Input file " << fileName << " could not be found!!" << endl;
	}
	
	int l=0;
  node_count = 0;
  elem_count = 0;
  
  dim = 3;
  
  bool start_node = false;
  bool start_elem = false;
  
  int line_start_node;
	int line_start_elem;
  
	if (found) {	
		while(getline(file, line)) {
      rawData.push_back(line /*+ "\n"*/); 
      //if (strcmp(str_inp1, str_inp2) == 0)
      //Increment nodes
      //if (strcmp(str_inp1, str_inp2) == 0)
        //or str.compare
      //cout << "Searching "<<line.substr(0,4)<<endl;
      if (line.substr(0,4) == string("GRID")){
        //cout << "Node found!"<<endl;
        node_count++;
        if (!start_node){
          start_node = true;
          line_start_node = l;
        }
      } else if (line.substr(0,5) == string("CTRIA") || line.substr(0,5) == string("CBEAM") || line.substr(0,4) == string("CBAR") ){
        if (!start_elem){
          start_elem = true;
					line_start_elem = l;
				}
        if (line.substr(0,5) == string("CBEAM") || line.substr(0,4) == string("CBAR"))
          dim = 2;
        //cout << "Element found!"<<endl;
        elem_count++;
      }
      l++;
    }
		file.close();
     
    cout << node_count <<" nodes and "<<elem_count<< " elements found."<<endl;

		// Strip all the inline or block comments (C++ style) from the rawData
		//stripComments(rawData);
		// Strip all the white spaces from the rawData
		//strip_white_spaces(rawData);
	}
	cout << "[I] "<<l << " lines readed ..." <<endl;
	line_count = l;
  
  //Allocating nodes 
  cout << "Allocating nodes"<<endl;
  node  	= new double 	[3 * node_count];
  nodeid  = new int 		[node_count];

	// NODAL FIELD DATA IS: GRID|ID|CP|X1|	
  int curr_line = line_start_node;
	l = curr_line;
	Vec3_t min( 1000., 1000., 1000.);
  Vec3_t max(-1000.,-1000.,-1000.);
	
	for (int n=0;n<node_count;n++){
    //cout << n+1; //DEBUG
		string temp = rawData[l].substr(FIELD_LENGTH,FIELD_LENGTH); //Second field, id
		nodeid[n] = atoi(temp.c_str());
		nodepos.insert(std::make_pair(atoi(temp.c_str()),n));
		//cout << "id: "<<nodeid[n]<<endl;
		for (int i = 0;i<3;i++) {
			int pos = 3*(FIELD_LENGTH)+ i*FIELD_LENGTH;
			//cout << "pos: "<<pos<<endl; 
			string temp = rawData[l].substr(pos,FIELD_LENGTH);
			
			//Search + or - symbol (SCIENTIFIC NOTATION)
			//BUT! If these signs (mainly the minus), are the FIRST character there is 
			//not the exponential sign
			int sign_pos = 0;
			if 			(temp.find("+")!=std::string::npos) {
				sign_pos = temp.find("+"); //Performance....
				temp.insert(sign_pos,"E");
			}
			else if (temp.find("-")!=std::string::npos) {
				sign_pos = temp.find("-");
				if (sign_pos!=0)
					temp.insert(sign_pos,"E");
				else { //Search for another "-" (NOW IT WILL BE SCIENTIFIC NOTATION)
					sign_pos = temp.find("-",1);
					if (sign_pos !=std::string::npos) {
						temp.insert(sign_pos,"E");
					}
				}
			}	
			
			double d = strtod(temp.c_str(),NULL);
			//cout << temp<<", conv: "<<d<<"sign pos" << sign_pos<<endl;
			//cout <<d<< " ";
			node[3*n+i] = d;
			if (d<min[i])
				min[i] = d;
			else if (d > max[i])
				max[i] = d;
		}
    //cout <<endl;
		l++;
  }
	
	cout << "Min values: "<< min <<endl;
	cout << "Max values: "<< max <<endl;	
  
  //IF FIXED FIELD
  cout << "Allocating Elements..."<<endl;
	// ASSUMING NODE IS FROM 1
  elcon = new int    [3 * elem_count];

	map<int, int>::iterator it;
  curr_line = line_start_elem;
	l = curr_line;
  for (int n=0;n<elem_count;n++){
    //cout << n+1<< " ";
		for (int en=0;en<dim;en++){
			int pos = 3*(FIELD_LENGTH)+ en*FIELD_LENGTH;
			string temp = rawData[l].substr(pos,FIELD_LENGTH); //Second field, id
			int d = atoi(temp.c_str());
			int nod = nodepos.find(d)->second;
			//cout << "node ind: "<<d<<"real node ind: "<<nod<<endl; 
			elcon[3*n+en] = nod;
			//cout << d<<" ";
		}
		//cout << endl;
		l++;
	}    
  cout << "Done."<<endl;
	
}

void NastranReader::WriteCSV(char const * FileKey)
{
	//type definition to shorten coding
	std::ostringstream oss;
	//Writing in a Log file
	String fn(FileKey);
	
	oss << "X, Y, Z"<<endl;;
	
	//#pragma omp parallel for schedule(static) num_threads(Nproc)
	// #ifdef __GNUC__
	// for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
	// #else
	for (int i=0; i<node_count; i++)//Like in Domain::Move
	//#endif
	{
		for (int j=0;j<3;j++){
			oss << node[3*i+j];
			if (j<2)	oss <<", ";
		}
		oss <<endl;
		
	}

	fn = FileKey;
	fn.append(".csv");	
	std::ofstream of(fn.CStr(), std::ios::out);
	of << oss.str();
	of.close();
}

void NastranReader::WriteVTK(char const * FileKey)
{

	string fileName(FileKey);
	fileName.append(".vtu");	
	ofstream file;
	file.open((fileName).c_str(),ios::out);
	file << "<?xml version=\"1.0\"?>" << endl;
	file << "<VTKFile type=\"UnstructuredGrid\">" << endl;
	file << "<UnstructuredGrid>" << endl;
	file << "<Piece NumberOfPoints=\"" << node_count << "\" NumberOfCells=\"" << elem_count << "\">" << endl;
	file << "<Points>" << endl;
	file << "<DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" >" << endl;
	for (int n=0;n<node_count;++n) {
		for (int i=0; i<3; ++i) 
			file<< setw(16) << setprecision(8) << scientific << node[3*n+i] << endl;
	}
	file << "</DataArray>" << endl;
	file << "</Points>" << endl;
	file << "<Cells>" << endl;

	file << "<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\" >" << endl;
	for (int c=0;c<elem_count;++c) {
		for (int n=0;n<3;++n) {
			file << elcon[3*c+n] << "\t";
		}
		file << endl;
	}

	file << "</DataArray>" << endl;
	file << "<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\" >" << endl;
	int offset=0;
	for (int c=0;c<elem_count;++c) {
		offset+=3;
		file << offset << endl;
	}
	file << "</DataArray>" << endl;

	file << "<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\" >" << endl;
	for (int c=0;c<elem_count;++c) {
		file << "5" << endl; // Tetra //Code is "5 " OR "VTK_TRIANGLE"
		// if (Cell(c).Num_Vertex()==4) file << "10" << endl; // Tetra
		// if (Cell(c).Num_Vertex()==8) file << "12" << endl; // Hexa
		// if (Cell(c).Num_Vertex()==6) file << "13" << endl; // Prism
		// if (Cell(c).Num_Vertex()==5) file << "14" << endl; // Pyramid (Wedge)
	}
	file << endl;
	file << "</DataArray>" << endl;;

	file << "</Cells>" << endl;

	file << "<PointData Scalars=\"scalars\" format=\"ascii\">" << endl;


	file << "</PointData>" << endl;

	file << "</Piece>" << endl;
	file << "</UnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;
	file.close();

	return;
}

NastranReader::~NastranReader(){
  
  delete node;
  delete elcon;  
}

};

#endif

