

#include <Omega_h_amr.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_map.hpp>
#ifdef CUDA_BUILD
#include <cuda_runtime.h>
#endif

void create_mesh(Omega_h::Mesh& mesh, 
#ifdef CUDA_BUILD
                 double* d_node_coords, int num_nodes, 
                 int* d_element_conn, int num_elements
#else
                 double* h_node_coords, int num_nodes, 
                 int* h_element_conn, int num_elements
#endif
                 ) 
{
#ifdef CUDA_BUILD
    // GPU Case: Use Omega_h::Write<> with device pointers
    Omega_h::Write<Omega_h::Real> device_coords(d_node_coords, num_nodes * 3);
    Omega_h::Write<Omega_h::LO> device_tets(d_element_conn, num_elements * 4);
#else
    // CPU Case: Copy raw pointer data to Omega_h::HostWrite<>
    Omega_h::HostWrite<Omega_h::Real> coords(num_nodes * 3);
    Omega_h::HostWrite<Omega_h::LO> tets(num_elements * 4);

    for (int i = 0; i < num_nodes * 3; ++i) coords[i] = h_node_coords[i];
    for (int i = 0; i < num_elements * 4; ++i) tets[i] = h_element_conn[i];

    // Convert HostWrite to Write<> for Omega_h
    Omega_h::Write<Omega_h::Real> device_coords(coords);
    Omega_h::Write<Omega_h::LO> device_tets(tets);
#endif

    // Build mesh (works on both CPU and GPU)
    mesh.build_from_elems(Omega_h::OMEGA_H_SIMPLEX, 3, device_tets, device_coords);
}

int main(int argc, char** argv) {
    Omega_h::Library lib;
    Omega_h::Mesh mesh(&lib);

    // Example data
    int num_nodes = 4;
    int num_elements = 1;

    double h_node_data[] = {
        0.0, 0.0, 0.0,  // v0
        1.0, 0.0, 0.0,  // v1
        0.0, 1.0, 0.0,  // v2
        0.0, 0.0, 1.0   // v3
    };

    int h_connectivity_data[] = {
        0, 1, 2, 3  // One tetrahedron
    };

#ifdef CUDA_BUILD
    // Allocate GPU memory
    double* d_node_data;
    int* d_connectivity_data;
    cudaMalloc((void**)&d_node_data, num_nodes * 3 * sizeof(double));
    cudaMalloc((void**)&d_connectivity_data, num_elements * 4 * sizeof(int));

    // Copy from host to GPU
    cudaMemcpy(d_node_data, h_node_data, num_nodes * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_connectivity_data, h_connectivity_data, num_elements * 4 * sizeof(int), cudaMemcpyHostToDevice);

    // Create mesh using GPU data
    create_mesh(mesh, d_node_data, num_nodes, d_connectivity_data, num_elements);

    // Free GPU memory
    cudaFree(d_node_data);
    cudaFree(d_connectivity_data);
#else
    // CPU case
    create_mesh(mesh, h_node_data, num_nodes, h_connectivity_data, num_elements);
#endif

    // Save mesh
    Omega_h::write_mesh("output.osh", &mesh);
    Omega_h::vtk_export_mesh("output.vtk", &mesh);

    return 0;
}
