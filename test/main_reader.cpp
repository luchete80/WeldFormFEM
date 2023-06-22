#include "NastranReader.h"

int main (){
  
  double **node; 
  int **elcon;
  ReadNastranTriMesh( "cylinder.nas", node, elcon);
}