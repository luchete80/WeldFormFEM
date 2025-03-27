#include <array>
#include <iostream>

///FOR THE ERROR  
//std::array<double, 2> target_node = {new_mesh.coords()[n][0], new_mesh.coords()[n][1]};
//double target_node[2] = {new_mesh.coords()[n][0], new_mesh.coords()[n][1]};


// Function to compute barycentric coordinates for a 3D tetrahedron
std::array<double, 4> barycentric_coordinates(const std::array<double, 3>& p,
                                              const std::array<double, 3>& p0,
                                              const std::array<double, 3>& p1,
                                              const std::array<double, 3>& p2,
                                              const std::array<double, 3>& p3) {
    // Compute volume of the tetrahedron
    std::array<double, 3> v0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    std::array<double, 3> v1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    std::array<double, 3> v2 = {p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]};
    
    double detT = v0[0] * (v1[1] * v2[2] - v1[2] * v2[1])
                - v0[1] * (v1[0] * v2[2] - v1[2] * v2[0])
                + v0[2] * (v1[0] * v2[1] - v1[1] * v2[0]);

    if (detT == 0.0) {
        std::cerr << "Degenerate tetrahedron encountered!" << std::endl;
        return {-1.0, -1.0, -1.0, -1.0}; // Invalid barycentric coordinates
    }

    // Compute barycentric coordinates
    std::array<double, 3> vp = {p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]};

    double lambda0 = ((vp[0] * (v1[1] * v2[2] - v1[2] * v2[1]))
                    - (vp[1] * (v1[0] * v2[2] - v1[2] * v2[0]))
                    + (vp[2] * (v1[0] * v2[1] - v1[1] * v2[0]))) / detT;

    double lambda1 = ((v0[0] * (vp[1] * v2[2] - vp[2] * v2[1]))
                    - (v0[1] * (vp[0] * v2[2] - vp[2] * v2[0]))
                    + (v0[2] * (vp[0] * v2[1] - vp[1] * v2[0]))) / detT;

    double lambda2 = ((v0[0] * (v1[1] * vp[2] - v1[2] * vp[1]))
                    - (v0[1] * (v1[0] * vp[2] - v1[2] * vp[0]))
                    + (v0[2] * (v1[0] * vp[1] - v1[1] * vp[0]))) / detT;

    double lambda3 = 1.0 - lambda0 - lambda1 - lambda2;

    return {lambda0, lambda1, lambda2, lambda3};
}

// Function to interpolate scalar values at the nodes of a tetrahedron
double interpolate_scalar(const std::array<double, 3>& p,
                          const std::array<double, 3>& p0, const std::array<double, 3>& p1, 
                          const std::array<double, 3>& p2, const std::array<double, 3>& p3,
                          double scalar0, double scalar1, double scalar2, double scalar3) {
    auto lambdas = barycentric_coordinates(p, p0, p1, p2, p3);
    return lambdas[0] * scalar0 + lambdas[1] * scalar1 + lambdas[2] * scalar2 + lambdas[3] * scalar3;
}

// Function to interpolate vector values at the nodes of a tetrahedron
std::array<double, 3> interpolate_vector(const std::array<double, 3>& p,
                                         const std::array<double, 3>& p0, const std::array<double, 3>& p1,
                                         const std::array<double, 3>& p2, const std::array<double, 3>& p3,
                                         std::array<double, 3> v0, std::array<double, 3> v1,
                                         std::array<double, 3> v2, std::array<double, 3> v3) {
    auto lambdas = barycentric_coordinates(p, p0, p1, p2, p3);
    return {
        lambdas[0] * v0[0] + lambdas[1] * v1[0] + lambdas[2] * v2[0] + lambdas[3] * v3[0],
        lambdas[0] * v0[1] + lambdas[1] * v1[1] + lambdas[2] * v2[1] + lambdas[3] * v3[1],
        lambdas[0] * v0[2] + lambdas[1] * v1[2] + lambdas[2] * v2[2] + lambdas[3] * v3[2]
    };
}

template <int dim>
void Map(Mesh& mesh, const Mesh& new_mesh) {
    // Loop over the target nodes in the new mesh
    for (int n = 0; n < new_mesh.nverts(); n++) {
        bool found = false;  // Flag to indicate whether the node is inside an element in the old mesh

        // Get coordinates for the node in the new mesh
        std::array<double, 2> target_node = {new_mesh.coords()[n][0], new_mesh.coords()[n][1]};

        // Loop over the elements in the old mesh (using *elnod to access connectivity and *node for coordinates)
        for (int i = 0; i < num_elements_in_old_mesh; i++) {
            // Connectivity for the tetrahedral element (assumed to have 4 nodes per element in the old mesh)
            int n0 = elnod[4*i];   // Node 0 in the element
            int n1 = elnod[4*i+1]; // Node 1 in the element
            int n2 = elnod[4*i+2]; // Node 2 in the element
            int n3 = elnod[4*i+3]; // Node 3 in the element
            
            // Coordinates of the element nodes (in the old mesh)
            std::array<double, 2> p0 = {node[2*n0], node[2*n0+1]};  // 2D coordinates of node 0
            std::array<double, 2> p1 = {node[2*n1], node[2*n1+1]};  // 2D coordinates of node 1
            std::array<double, 2> p2 = {node[2*n2], node[2*n2+1]};  // 2D coordinates of node 2
            std::array<double, 2> p3 = {node[2*n3], node[2*n3+1]};  // 2D coordinates of node 3

            // Check if the target node is inside the element (using barycentric coordinates for 2D triangle)
            std::array<double, 3> lambdas = barycentric_coordinates(target_node, p0, p1, p2);

            // If the target node is inside the triangle (element)
            if (lambdas[0] >= -5.0e-2 && lambdas[1] >= -5.0e-2 && lambdas[2] >= -5.0e-2) {
                // Interpolate scalar values at the element's nodes
                double scalar0 = 0.0;  // Replace with actual scalar value at node 0 in the old mesh
                double scalar1 = 0.0;  // Replace with actual scalar value at node 1 in the old mesh
                double scalar2 = 0.0;  // Replace with actual scalar value at node 2 in the old mesh
                double interpolated_scalar = interpolate_scalar(target_node, p0, p1, p2, scalar0, scalar1, scalar2);

                // Interpolate vector values for displacement (if needed)
                std::array<double, 3> disp0 = {0.0, 0.0, 0.0};  // Replace with displacement at node 0
                std::array<double, 3> disp1 = {0.0, 0.0, 0.0};  // Replace with displacement at node 1
                std::array<double, 3> disp2 = {0.0, 0.0, 0.0};  // Replace with displacement at node 2
                std::array<double, 3> interpolated_disp = interpolate_vector(target_node, p0, p1, p2, disp0, disp1, disp2);

                // Optionally, interpolate other scalar/vector fields for the new mesh node here

                std::cout << "Node " << n << " is inside element " << i << " of the old mesh." << std::endl;
                std::cout << "Interpolated scalar: " << interpolated_scalar << std::endl;
                std::cout << "Interpolated displacement: (" << interpolated_disp[0] << ", " << interpolated_disp[1] << ", " << interpolated_disp[2] << ")\n";

                found = true;
                break;  // Exit the element loop once the element is found
            }
        }

        if (!found) {
            std::cout << "Node " << n << " is not inside any element of the old mesh." << std::endl;
        }
    }
}