/*************************************************************************/
/*  test_omega_h_adapt.cpp                                       */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                              */
/*  https://www.opensourcemech.com                                                                */
/*                                                                       */
/*  Copyright (c) 2025-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/


#include <Omega_h.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_vtk.hpp>

#ifdef CUDA_BUILD
#include <cuda_runtime.h>
#endif

// Function to evaluate if refinement is necessary based on the scalar field (e.g., plastic strain)
bool needs_refinement(double scalar_value, double threshold) {
    return scalar_value > threshold;
}

// Refinement function (simplified version for illustration)
void refine_mesh(Omega_h::Mesh& mesh, double threshold, 
                 Omega_h::Write<Omega_h::Real>& plastic_strain_field) {
    // Loop through elements or nodes to check the plastic strain and refine mesh
    // Assuming the scalar field is stored at the element level or node level

    // For this example, we're checking at the element level
    for (int i = 0; i < mesh.nents(0); ++i) {
        Omega_h::LOs elem_nodes = mesh.conn(0, i);  // Get nodes of the element
        double max_strain = 0.0;

        // Evaluate the maximum strain in the element (simplified)
        for (auto node : elem_nodes) {
            double strain = plastic_strain_field[node];
            max_strain = std::max(max_strain, strain);
        }

        // If the strain exceeds the threshold, refine the element
        if (needs_refinement(max_strain, threshold)) {
            // Implement refinement (this can be splitting the element, adding new nodes, etc.)
            // Here we will just print it for demonstration, but in practice,
            // you'd refine the element (e.g., split or add new nodes)
            std::cout << "Element " << i << " needs refinement due to strain: " << max_strain << "\n";
            // Example: mesh.refine_element(i); // Hypothetical function for refinement
        }
    }
}

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

    // Example plastic strain field (this would typically come from a simulation)
    double plastic_strain_data[] = {0.1, 0.2, 0.6, 0.8};  // Example plastic strain at each node

    // Refinement threshold for plastic strain
    double strain_threshold = 0.5;

#ifdef CUDA_BUILD
    // Allocate GPU memory
    double* d_node_data;
    int* d_connectivity_data;
    double* d_plastic_strain_data;
    cudaMalloc((void**)&d_node_data, num_nodes * 3 * sizeof(double));
    cudaMalloc((void**)&d_connectivity_data, num_elements * 4 * sizeof(int));
    cudaMalloc((void**)&d_plastic_strain_data, num_nodes * sizeof(double));

    // Copy data from host to GPU
    cudaMemcpy(d_node_data, h_node_data, num_nodes * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_connectivity_data, h_connectivity_data, num_elements * 4 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_plastic_strain_data, plastic_strain_data, num_nodes * sizeof(double), cudaMemcpyHostToDevice);

    // Create mesh using GPU data
    create_mesh(mesh, d_node_data, num_nodes, d_connectivity_data, num_elements);

    // Create plastic strain field (using GPU data)
    Omega_h::Write<Omega_h::Real> plastic_strain_field(d_plastic_strain_data, num_nodes);

    // Refine mesh based on plastic strain
    refine_mesh(mesh, strain_threshold, plastic_strain_field);

    // Free GPU memory
    cudaFree(d_node_data);
    cudaFree(d_connectivity_data);
    cudaFree(d_plastic_strain_data);
#else
    // Create plastic strain field (using CPU data)
    Omega_h::Write<Omega_h::Real> plastic_strain_field(plastic_strain_data, num_nodes);

    // Create mesh using CPU data
    create_mesh(mesh, h_node_data, num_nodes, h_connectivity_data, num_elements);

    // Refine mesh based on plastic strain
    refine_mesh(mesh, strain_threshold, plastic_strain_field);
#endif

    // Export to VTK
    Omega_h::vtk_export_mesh("output_refined.vtk", &mesh);

    return 0;
}
