#include "Domain.h"

#include <iostream>
using namespace std;

int main(){
	
	Domain dom(2);

	Vec3D pos(0.,0.,0.);
	Vec3D L(10.0,1.0,0.);
	double r = 0.1;
	
	dom.AddBoxLength(pos,L,r);
	
	cout << "Program ended."<<endl;
  
  return 0;
  
}