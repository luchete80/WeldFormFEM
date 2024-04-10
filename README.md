# EuLaFEM
Eulerian Lagrange FEM



CUDACXX=/usr/local/cuda-12.3/bin/nvcc cmake ../WeldFormFEM -DBUILD_GPU=ON

To update libraries (LSDynaReader and Math)

git submodule update --init --recursive
