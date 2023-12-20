#include "NastranReader.h"

int main (){
  
  double *node = NULL;
  //*node  = new double[1]; 
  int *elcon =NULL;
  int nodecount;
  ReadNastranTriMesh( "cylinder.nas", &node, &elcon,&nodecount); //This works calling by the address
}