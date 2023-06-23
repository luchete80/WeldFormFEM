#include "NastranReader.h"

int main (){
  
  double *node = NULL;
  //*node  = new double[1]; 
  int *elcon =NULL;
  ReadNastranTriMesh( "cylinder.nas", &node, &elcon); //This works calling by the address
}