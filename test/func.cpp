#include <iostream>
#include "addition.h"
using namespace std;

#ifdef __cplusplus
  extern"C" {
#endif

void *__gxx_personality_v0;

int func(double a[]) {
   int i;
   for(i=0; i<4; i++) {
       cout << a[i] << endl;
   }
   int z;
   z = addition (5,3);
   cout << z << endl;
   return 0;
}

#ifdef __cplusplus
  }
#endif