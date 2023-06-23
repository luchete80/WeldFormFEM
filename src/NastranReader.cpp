#include "NastranReader.h"

extern "C" void ReadNastranTriMesh( char* fName, int *nodecount, double **node, int **elcon)
{
  
  std::map <int,int> nodepos;	//id to position
	// std::vector <std::string> rawData;
	// int line_count;
  
	string fileName = fName;
  string line;
  //rawData="";
	fstream file;
    bool found=false;
	//MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	cout << "[I] Reading ... "<<fileName.c_str()<<endl;
	file.open(fileName.c_str());
	if (file.is_open()) {
		cout << "[I] Found input file " << fileName << endl;
		found=true;
	} else {
		cerr << "[E] Input file " << fileName << " could not be found!!" << endl;
	}
	
	int l=0;
  int nodxelem;
  int node_count = 0;
  int elem_count = 0;
	std::vector <std::string> rawData;
  int *nodeid, line_count;
  
  bool issurf = false;
  int dim = 3;
  
  //  Still pending with tri and tetra nonsurface element 
  if (dim == 3){
    if (issurf) nodxelem = 3;
    else        nodxelem = 8;
  } else if (dim == 2){
    if (issurf) nodxelem = 2;
    else        nodxelem = 4;    
  }
  
  cout << "Dimension: "<<dim <<", Is surface Mesh: "<<issurf<<endl; 
  cout << "Nodes per element: " <<nodxelem<<endl;
  
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
      } else if ( line.substr(0,5) == string("CTRIA") || line.substr(0,5) == string("CBEAM") || line.substr(0,4) == string("CBAR")
               || line.substr(0,5) == string("CQUAD") || line.substr(0,5) == string("CHEXA")){ ////// SURFMESH
        if (!start_elem){
          start_elem = true;
					line_start_elem = l;
				}
        // if (line.substr(0,5) == string("CBEAM") || line.substr(0,4) == string("CBAR"))
          // dim = 2;
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
  *nodecount = node_count;
	cout << "[I] "<<l << " lines readed ..." <<endl;
	line_count = l;
  
  //Allocating nodes 
  cout << "Allocating nodes"<<endl;
  *node  	= new double 	[nodxelem * node_count];
  nodeid  = new int 		[node_count];

	// NODAL FIELD DATA IS: GRID|ID|CP|X1|	
  int curr_line = line_start_node;
	l = curr_line;
	//Vec3_t min( 1000., 1000., 1000.);
  //Vec3_t max(-1000.,-1000.,-1000.);
	//cout << "Reading Nodes  ... "<<endl;
	for (int n=0;n<node_count;n++){
    //cout << n+1; //DEBUG
		string temp = rawData[l].substr(FIELD_LENGTH,FIELD_LENGTH); //Second field, id
    //cout << "temp "<<temp<<endl;
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
			(*node)[nodxelem*n+i] = d;
			// if (d<min[i])
				// min[i] = d;
			// else if (d > max[i])
				// max[i] = d;
		}
		l++;
  }
	
	// cout << "Min values: "<< min <<endl;
	// cout << "Max values: " << max <<endl;	
  
  //IF FIXED FIELD
  cout << "Allocating Elements..."<<endl;
	// ASSUMING NODE IS FROM 1
  *elcon = new int    [nodxelem * elem_count];

	map<int, int>::iterator it;
  curr_line = line_start_elem;
	l = curr_line;
  int fieldnum[] ={6,6,8}; //per line
  for (int n=0;n<elem_count;n++){
    //cout << n+1<< " ";
    if (issurf){
      
      for (int en=0;en<dim;en++){
        int pos = nodxelem*(FIELD_LENGTH)+ en*FIELD_LENGTH;
        string temp = rawData[l].substr(pos,FIELD_LENGTH); //Second field, id
        int d = atoi(temp.c_str());
        int nod = nodepos.find(d)->second;
        //cout << "node ind: "<<d<<"real node ind: "<<nod<<endl; 
        (*elcon)[nodxelem*n+en] = nod;
        //cout << d<<" ";
      }
		} else  {
      for (int lin=0;lin<3;lin++){
        
      }
    }
    //cout << endl;
		l++;
	}    
  cout << "Done."<<endl;
	
}


