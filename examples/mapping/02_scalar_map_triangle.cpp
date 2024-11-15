#include <triangle.h>
#include <vector>
#include <iostream>
#include <array>

// Function to compute barycentric coordinates
std::array<double, 3> barycentric_coordinates(const std::array<double, 2>& p,
                                              const std::array<double, 2>& p0,
                                              const std::array<double, 2>& p1,
                                              const std::array<double, 2>& p2) {
    double denominator = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
    double lambda1 = ((p1[0] - p[0]) * (p2[1] - p[1]) - (p2[0] - p[0]) * (p1[1] - p[1])) / denominator;
    double lambda2 = ((p2[0] - p[0]) * (p0[1] - p[1]) - (p0[0] - p[0]) * (p2[1] - p[1])) / denominator;
    double lambda3 = 1.0 - lambda1 - lambda2;
    return {lambda1, lambda2, lambda3};
}

// Function to interpolate scalar values
double interpolate_scalar(const std::array<double, 2>& p,
                          const std::array<double, 2>& p0, const std::array<double, 2>& p1, const std::array<double, 2>& p2,
                          double scalar0, double scalar1, double scalar2) {
    auto lambdas = barycentric_coordinates(p, p0, p1, p2);
    return lambdas[0] * scalar0 + lambdas[1] * scalar1 + lambdas[2] * scalar2;
}

int main() {
    // Define the source mesh
    std::vector<std::array<double, 2>> source_nodes = {
        {0, 0}, {1, 0}, {0, 1}, {1, 1}
    };
    std::vector<std::array<int, 3>> source_triangles = {
        {0, 1, 2}, {1, 3, 2}
    };
    std::vector<double> source_scalars = {1.0, 2.0, 3.0, 4.0};

    // Define the target points
    std::vector<std::array<double, 2>> target_points = {
        {0.25, 0.25}, {0.75, 0.25}, {0.25, 0.75}, {0.5, 0.5}
    };

    // Prepare Triangle data structures
    struct triangulateio in, out;
    memset(&in, 0, sizeof(in));
    memset(&out, 0, sizeof(out));

    // Load source nodes into Triangle
    in.numberofpoints = source_nodes.size();
    in.pointlist = (double*)malloc(in.numberofpoints * 2 * sizeof(double));
    for (size_t i = 0; i < source_nodes.size(); ++i) {
        in.pointlist[2 * i] = source_nodes[i][0];
        in.pointlist[2 * i + 1] = source_nodes[i][1];
    }

    // Load source triangles into Triangle
    in.numberoftriangles = source_triangles.size();
    in.trianglelist = (int*)malloc(in.numberoftriangles * 3 * sizeof(int));
    for (size_t i = 0; i < source_triangles.size(); ++i) {
        in.trianglelist[3 * i] = source_triangles[i][0];
        in.trianglelist[3 * i + 1] = source_triangles[i][1];
        in.trianglelist[3 * i + 2] = source_triangles[i][2];
    }

    // Perform triangulation
    triangulate("z", &in, &out, nullptr);

    // Interpolate scalar values for each target point
    for (const auto& target : target_points) {
        bool found = false;
        for (int i = 0; i < out.numberoftriangles; ++i) {
            int v0 = out.trianglelist[3 * i];
            int v1 = out.trianglelist[3 * i + 1];
            int v2 = out.trianglelist[3 * i + 2];

            std::array<double, 2> p0 = {out.pointlist[2 * v0], out.pointlist[2 * v0 + 1]};
            std::array<double, 2> p1 = {out.pointlist[2 * v1], out.pointlist[2 * v1 + 1]};
            std::array<double, 2> p2 = {out.pointlist[2 * v2], out.pointlist[2 * v2 + 1]};

            // Check if the point is inside the triangle
            auto lambdas = barycentric_coordinates(target, p0, p1, p2);
            if (lambdas[0] >= 0 && lambdas[1] >= 0 && lambdas[2] >= 0) {
                double scalar0 = source_scalars[v0];
                double scalar1 = source_scalars[v1];
                double scalar2 = source_scalars[v2];
                double interpolated_value = interpolate_scalar(target, p0, p1, p2, scalar0, scalar1, scalar2);
                std::cout << "Interpolated scalar at (" << target[0] << ", " << target[1] << ") = " 
                          << interpolated_value << "\n";
                found = true;
                break;
            }
        }
        if (!found) {
            std::cout << "Point (" << target[0] << ", " << target[1] << ") is outside the mesh.\n";
        }
    }

    // Clean up
    free(in.pointlist);
    free(in.trianglelist);
    free(out.pointlist);
    free(out.trianglelist);

    return 0;
}