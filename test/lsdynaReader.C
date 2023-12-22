#include "lsdynaReader.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdio.h>

using namespace std;

bool isComment(const string &line){
  
  bool ret = false;
  if (line.find('$') == 0)
    ret = true;
  return ret;
}
// str.erase(
  // std::remove_if(str.begin(), str.end(), [](char chr){ return chr == '&' || chr == ' ';}),
  // str.end());

int findNextCommandLine(const int &curr_line, std::vector <std::string> line) {
  string first;
  int i = curr_line;
  bool end = false;
  while (!end) {
    i++;
    if (i==line.size()) end = true; 
      first = line[i].substr(0,1);
      //cout << "first "<<first <<", "<<line[i]<<endl;
      //Can a line start with spaces??
      //if (first.find("$")!=string::npos)
      if (first=="*"){
        //cout << "FOUND * at line "<<i<<endl;
        end = true;
      }
  }// while end
  
        //cout << "* NOT FOUND "<<i<<endl;
  return i;
}

void removeComments(std::vector <string> m_line){
  std::vector<std::string>::iterator it;
  int i=0;
  for (it = m_line.begin();it!=m_line.end();){
    //cout << *it<<endl;
    string first = it->substr(0,1);
    //cout << first <<endl;
    if (first=="$"){
      m_line.erase(it);
      //m_line_count--;
    } else {
      it ++;
      //cout << "NOT $ FOUND "<<first<<", LIN: "<<*it<<endl;
    }
    i++;
  }
}

void lsdynaReader::removeComments(){
  std::vector<std::string>::iterator it;
  int i=0;
  for (it = m_line.begin();it!=m_line.end();){
    //cout << *it<<endl;
    string first = it->substr(0,1);
    //cout << first <<endl;
    if (first=="$"){
      m_line.erase(it);
      m_line_count--;
    } else {
      it ++;
      //cout << "NOT $ FOUND "<<first<<", LIN: "<<*it<<endl;
    }
    i++;
  }
}

string removeSpaces(string str) { 
    str.erase(remove(str.begin(), str.end(), ' '), str.end()); 
    return str; 
}

double readDoubleField(string &str, const int &pos, const int &length) {
  double ret = 0.0;
    std::string f_str = str.substr(pos,length);
    f_str = removeSpaces(f_str);

  if (str.size()<pos+length) {
    cout << "ERROR TRYING TO READ: "<<str<<", subs: " <<f_str<<endl;
    return ret;
  }
  // return ret;
  ret = strtod(f_str.c_str(),NULL);
  // cout << "readed: "<<ret;
  return ret;
}

int readIntField(string &str, const int &pos, const int &length) {
  int ret = 0.0;
    std::string f_str = str.substr(pos,length);
    f_str = removeSpaces(f_str);

  if (str.size()<pos+length) {
    cout << "ERROR TRYING TO READ: "<<str<<", subs: " <<f_str<<endl;
    return ret;
  }
  // return ret;
  ret = stoi(f_str.c_str());
  //cout << "readed string "<<f_str<<endl;
  // cout << "readed: "<<ret;
  return ret;
}

// bool findSection(char **, string str, int * ini_pos, int *end_pos){
  
  // bool end = false;
  // bool found = false;
  // int i = 0;
  // cout << "Reading Elements"<<endl;
  // while (!end){

      // if (m_line[i].find(str) != std::string::npos){
        // found = true;
        // cout << "Found section" << str << " at line "<<i<<endl;
        // *ini_pos = i+1;
        // m_elem_count = findNextCommandLine(i,m_line) - i;
        // *end_pos = *ini_pos + m_elem_count -1 ;
        // cout << "Section length: "<<m_elem_count<<endl;
        // end = true;
      // }
    // if (i==m_line_count) {
      // end = true;
      // cout << "ELEMENT not defined "<<endl;
    // }
    // i++;
  // } 
  // return found;
// }

bool lsdynaReader::findSection(string str, int * ini_pos, int *end_pos){
  
  bool end = false;
  bool found = false;
  int i = 0;
  cout << "Reading Elements"<<endl;
  while (!end){

      if (m_line[i].find(str) != std::string::npos){
        found = true;
        cout << "Found section" << str << " at line "<<i<<endl;
        *ini_pos = i+1;
        m_elem_count = findNextCommandLine(i,m_line) - i;
        *end_pos = *ini_pos + m_elem_count -1 ;
        cout << "Section length: "<<m_elem_count<<endl;
        end = true;
      }
    if (i==m_line_count) {
      end = true;
      cout << "ELEMENT not defined "<<endl;
    }
    i++;
  } 
  return found;
}

void lsdynaReader::readNodes() {
  bool end = false;
  int ini_pos, end_pos;
  int i = 0;
  findSection ("*NODE", &ini_pos, &end_pos);


  for (i=ini_pos;i<end_pos;i++){
    int id;
    id = readDoubleField(m_line[i], 0, 8);
    ls_node nod;
    nod.m_id = id;
    for (int d=0;d<3;d++)
      nod.m_x[d] = readDoubleField(m_line[i], 8+16*d, 16);
      // cout << "Node "<<id <<"XYZ: "<<nod.m_x[0]<<", "<<nod.m_x[1]<<", "<<nod.m_x[2]<<endl; 
    
    
  }

}  //line

void lsdynaReader::readElementSolid() {
  bool end = false;
  int ini_pos, end_pos;
  int i = 0;
  findSection ("*ELEMENT_SOLID", &ini_pos, &end_pos);
  
  for (i=ini_pos;i<end_pos;i++){
    int id, pid;
    ls_element ls_el;
    int nodecount;
    if (m_line[i].size()>16) nodecount = (int)((m_line[i].size()-16)/8);
    // cout << "Elem node count "<<nodecount <<endl;
    ls_el.id  = readIntField(m_line[i], 0, 8);
    ls_el.pid = readIntField(m_line[i], 1, 8);
    ls_node nod;
    nod.m_id = id;
    for (int d=0;d<nodecount;d++)
      ls_el.node.push_back(readIntField(m_line[i], 16+8*d, 8));
      // cout << "Node "<<id <<"XYZ: "<<nod.m_x[0]<<", "<<nod.m_x[1]<<", "<<nod.m_x[2]<<endl; 
      m_elem.push_back(ls_el);
    
  }

}  //line

/********************************************/
/* 1 NODE ID AND N DOF AND */
void lsdynaReader::readSPCNodes(){
  
    bool end = false;
  int ini_pos, end_pos;
  int i = 0;
  cout << "Searching Boundary SPC"<<endl;
  findSection ("*BOUNDARY_SPC_NODE", &ini_pos, &end_pos);
  
  int k = 0;
  for (i=ini_pos;i<end_pos;i++){
    int id, cid;
    ls_spc_node ls_nspc;
    // cout << "Elem node count "<<nodecount <<endl;
    ls_nspc.m_node_id  = readIntField(m_line[i], 0, 10);
    //cout << "m_node "<<ls_nspc.m_node_id<<endl;
    // ls_el.cid = readIntField(m_line[i], 1, 8);
    
    for (int d=0;d<6;d++)
      ls_nspc.m_fix_dof[d] = readIntField(m_line[i], 20+10*d, 10);
      // cout << "Node "<<id <<"XYZ: "<<nod.m_x[0]<<", "<<nod.m_x[1]<<", "<<nod.m_x[2]<<endl; 
    //cout << "spc: "<<k<<endl; k++;
    m_spc_nod.push_back(ls_nspc);    
  }
  
}

// $#     nid       dof       vad      lcid        sf       vid     death     birth
bool readBPMNodes() {
  
}

lsdynaReader::lsdynaReader(const char *fname){
  string line;
  m_line_count = 0;
  int start, end;
  char dl = ' ';
  ifstream file(fname);
  if (file.is_open()) {
    while (getline(file, line)) {
      m_line.push_back(line);
      m_line_count++;
    }
  }
        file.close();
  cout << "Line count: "<< m_line_count << endl;
  
  removeComments();
  cout << "Line count w/o comments: "<< m_line_count << endl;
  
  readNodes();
  readElementSolid();
  readSPCNodes();
  //CHECK FOR
// *BOUNDARY_SPC_SET
// $#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz
         // 1         0         1         1         1         1         1         1
// *SET_NODE_LIST_TITLE
// *BOUNDARY_PRESCRIBED_MOTION_NODE

}


//// NON CLASS /////
void readNodes(char *fName, double **nodes, int *node_count){
  string line;
  int m_line_count = 0;
  std::vector <string> m_line;
  int start, end;
  cout << "Opening:  "<<fName<< " ."<<endl;
  char dl = ' ';
  ifstream file(fName);
  if (file.is_open()) {
    while (getline(file, line)) {
      m_line.push_back(line);
      m_line_count++;
    }
  }
  cout << "lines: "<<m_line_count<<endl;
  *node_count = 100;
  
  // *nodes = (double *) malloc(3*reader.m_node_count*sizeof(double));
  *nodes = (double *) malloc(3*(*node_count)*sizeof(double));
  for (int n=0;n<*node_count;n++)
    for (int d=0;d<3;d++) {
    //(*nodes)[3*n+d] = reader.m_node[n].m_x[d];
    (*nodes)[3*n+d] = 1.;
  }
}

void LSDYNA_getLines(char* fname, char ***lines) {
  string line;
  std::vector <string> m_line;
  int m_line_count = 0;
  int start, end;
  char dl = ' ';
  ifstream file(fname);
  if (file.is_open()) {
    while (getline(file, line)) {
      m_line.push_back(line);
      m_line_count++;
    }
  }
  removeComments(m_line);
  
  cout << "Line size: "<< m_line_count<<endl;
  
  *lines = (char **) malloc(m_line_count * sizeof(char *));
  for (int l=0;l<m_line_count;l++){
    (*lines)[l] = (char *) malloc (100 * sizeof(char));
    //(*lines)[l] = (char *) m_line[l].c_str(); 
    for (int i=0;i<m_line[l].size();i++)
      (*lines)[l][i] = m_line[l].c_str()[i];
  }
  
}