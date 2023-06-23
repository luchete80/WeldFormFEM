#include "NastranReader.h"

//TEST
extern "C" void c_func(int **x, int s) {
  *x = (int *) malloc(s*sizeof(int));
  *x[0] = 100; (*x)[1] = 30; //atention, is not *x[] = 
}

//#ifdef _COMP_FORTRAN_ 
//extern "C" void ReadNastranTriMesh( char* fName, double **node, int **elcon)
//#else
  void ReadNastranTriMesh( char* fName, double **node, int **elcon)
{
  
 std::ofstream of("output.csv", std::ios::out);
std::ostringstream oss; 
 bool issurf = false;
  int node_count, elem_count, line_count;
  
  std::map <int,int> nodepos;	//id to position
  int dim;
	std::vector <std::string> rawData;
  int *nodeid;
	// int line_count;
  
	std::string fileName = fName;
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
  int nodxelem;

  if (dim == 3){
    if (issurf) nodxelem = 3;
    else        nodxelem = 8;
  } else if (dim == 2){
    if (issurf) nodxelem = 2;
    else        nodxelem = 4;    
  }
  
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
      } else if (line.substr(0,5) == string("CTRIA") || line.substr(0,5) == string("CBEAM") || line.substr(0,4) == string("CBAR") 
                     || line.substr(0,5) == string("CQUAD") || line.substr(0,5) == string("CHEXA")){
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
  //*node  	= new double 	[3 * node_count];
  *node = (double *) malloc(3*node_count*sizeof(double));
  nodeid  = new int 		[node_count];

	// NODAL FIELD DATA IS: GRID|ID|CP|X1|	
  int curr_line = line_start_node;
	l = curr_line;
  // double min[] = { 1000., 1000., 1000.};
  // double max[] = {-1000.,-1000.,-1000.};
	cout << "Reading nodes "<<endl;
	for (int n=0;n<node_count;n++){
    //cout << n+1; //DEBUG
		string temp = rawData[l].substr(FIELD_LENGTH,FIELD_LENGTH); //Second field, id
		//nodeid[n] = atoi(temp.c_str());
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
      
      of <<d;

      if (i<2) of<<", "; 
			(*node)[3*n+i] = d;
			// if (d<min[i])
				// min[i] = d;
			// else if (d > max[i])
				// max[i] = d;
		}//for node component i
    of<<endl;
		l++;
  }
	
	// cout << "Min values: "<< min <<endl;
	// cout << "Max values: "<< max <<endl;	
  
  //IF FIXED FIELD
  cout << "Allocating Elements..."<<endl;
	// ASSUMING NODE IS FROM 1
  *elcon = new int    [nodxelem * elem_count];

	map<int, int>::iterator it;
  curr_line = line_start_elem;
	l = curr_line;
  int fieldnum[] ={6,2}; //per line, in case of 20 nodes is : {6,6,8}
  string field;
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
      if (dim ==3){
      int init_pos;
      for (int lin=0;lin<2;lin++){
        if (lin==0) init_pos = 3 *  FIELD_LENGTH;
        else        init_pos = FIELD_LENGTH;
        for (int c=0;c<fieldnum[lin];c++){
          string temp = rawData[l+lin].substr(init_pos,FIELD_LENGTH);
          //cout << "line "<<temp;
          int d = atoi(temp.c_str());
          int nod = nodepos.find(d)->second;
          //cout << temp<<", number "<<d<<"; ";
          of << d;
          if (!(lin==1 && c==(fieldnum[lin]-1)))of << ", ";
          init_pos += FIELD_LENGTH;
        }
      }
      of <<endl;
      //cout << endl;
      }///dim=3
    } //NO IS SURF
    //cout << endl;
		if (dim==3 && !issurf) l+=2;
	} ///Eleem    
	of << oss.str();
	of.close();
  cout << "Done."<<endl;
	
} //function

