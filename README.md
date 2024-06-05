# EuLaFEM
Eulerian Lagrange FEM



CUDACXX=/usr/local/cuda-12.3/bin/nvcc cmake ../WeldFormFEM -DBUILD_GPU=ON

To update libraries (LSDynaReader and Math)

git submodule update --init --recursive



Link here #https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/ to see different architectures. 
