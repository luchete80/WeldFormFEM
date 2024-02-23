https://stackoverflow.com/questions/35603035/creating-an-object-in-device-code

#include <stdio.h>
#define N 4


class Polygon {
  protected:
    int width, height;
  public:
  __host__ __device__  void set_values (int a, int b)
      { width=a; height=b; }
  __host__ __device__  virtual int area ()
      { return 0; }
};

class Rectangle: public Polygon {
  public:
  __host__ __device__  int area ()
      { return width * height; }
};

class Triangle: public Polygon {
  public:
  __host__ __device__   int area ()
      { return (width * height / 2); }
};

__global__ void setup_f(Polygon ** d_polys) {
  int idx = threadIdx.x+blockDim.x*blockIdx.x;
  if (idx < N) {
    if (idx%2)
      d_polys[idx] = new Rectangle();
    else
      d_polys[idx] = new Triangle();
    d_polys[idx]->set_values(5,12);
}};

__global__ void area_f(Polygon ** d_polys) {
  int idx = threadIdx.x+blockDim.x*blockIdx.x;
  if (idx < N){
    printf("area of object %d = %d\n", idx, d_polys[idx]->area());
}};


int main () {

  Polygon **devPolys;
  cudaMalloc(&devPolys,N*sizeof(Polygon *));
  setup_f<<<1,N>>>(devPolys);
  area_f<<<1,N>>>(devPolys);
  cudaDeviceSynchronize();
}

$ nvcc -o t1086 t1086.cu
$ cuda-memcheck ./t1086
========= CUDA-MEMCHECK
area of object 0 = 30
area of object 1 = 60
area of object 2 = 30
area of object 3 = 60
========= ERROR SUMMARY: 0 errors
$