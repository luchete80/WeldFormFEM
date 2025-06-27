#include "Matrix.h"

int main(){
  Matrix A(3,3);
  A.Set(0,0, 3); A.Set(0,1, 2); A.Set(0,2, -1);
  A.Set(1,0, 2); A.Set(1,1, -2); A.Set(1,2, 4);
  A.Set(2,0, -1); A.Set(2,1, 0.5); A.Set(2,2, -1);

  Matrix b(3,1);
  b.Set(0,0, 1);
  b.Set(1,0, -2);
  b.Set(2,0, 0);


  // solution is [1,-2,-2]

  // solution is [1,-2,-2]

  //~ Matrix A(2,2); A.Set(0,0,2); A.Set(0,1,1); A.Set(1,0,1); A.Set(1,1,3);
  //~ Matrix b(2,1); b.Set(0,0,1); b.Set(1,0,2);
  
printf("A INVERSE\n");
  A.Inv().Print();
  printf("-------\n");
  // Suppose you invert A and multiply
  //Matrix x = A.Inv() * b;
  Matrix x(3,1);
  x = SolveLinearSystem(A, b);
  x.Print(); // Should show x ≈ [0; ~0.67]

}
