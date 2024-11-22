#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <vector>
#include <iostream>
#include <array>

// Define CGAL Kernel and Delaunay triangulation
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
typedef K::Point_2 Point;

// Function to compute barycentric coordinates
std::array<double, 3> barycentric_coordinates(const Point& p, const Point& p0, const Point& p1, const Point& p2) {
    double denominator = (p1.x() - p0.x()) * (p2.y() - p0.y()) - (p2.x() - p0.x()) * (p1.y() - p0.y());
    double lambda1 = ((p1.x() - p.x()) * (p2.y() - p.y()) - (p2.x() - p.x()) * (p1.y() - p.y())) / denominator;
    double lambda2 = ((p2.x() - p.x()) * (p0.y() - p.y()) - (p0.x() - p.x()) * (p2.y() - p.y())) / denominator;
    double lambda3 = 1.0 - lambda1 - lambda2;
    return {lambda1, lambda2, lambda3};
}

// Function to interpolate scalar value
double interpolate_scalar(const Point& p, const Point& p0, const Point& p1, const Point& p2, 
                          double scalar0, double scalar1, double scalar2) {
    auto lambdas = barycentric_coordinates(p, p0, p1, p2);
    return lambdas[0] * scalar0 + lambdas[1] * scalar1 + lambdas[2] * scalar2;
}

int main() {
    // Source mesh nodes and scalar values
    std::vector<Point> source_nodes = { Point(0, 0), Point(1, 0), Point(0, 1), Point(1, 1) };
    std::vector<std::array<int, 3>> source_triangles = { {0, 1, 2}, {1, 3, 2} };
    std::vector<double> source_scalars = { 1.0, 2.0, 3.0, 4.0 };

    // Target points
    std::vector<Point> target_points = { Point(0.25, 0.25), Point(0.75, 0.25), Point(0.25, 0.75), Point(0.5, 0.5) };

    // Perform Delaunay triangulation for spatial queries
    Delaunay triangulation;
    triangulation.insert(source_nodes.begin(), source_nodes.end());

    // Interpolate scalar values at target points
    for (const auto& target : target_points) {
        // Find the enclosing triangle
        auto face = triangulation.locate(target);
        if (triangulation.is_infinite(face)) {
            std::cout << "Point (" << target << ") is outside the mesh.\n";
            continue;
        }

        // Get triangle vertices
        auto vertex0 = face->vertex(0)->point();
        auto vertex1 = face->vertex(1)->point();
        auto vertex2 = face->vertex(2)->point();

        // Get scalar values for the triangle vertices
        int index0 = std::find(source_nodes.begin(), source_nodes.end(), vertex0) - source_nodes.begin();
        int index1 = std::find(source_nodes.begin(), source_nodes.end(), vertex1) - source_nodes.begin();
        int index2 = std::find(source_nodes.begin(), source_nodes.end(), vertex2) - source_nodes.begin();

        double scalar0 = source_scalars[index0];
        double scalar1 = source_scalars[index1];
        double scalar2 = source_scalars[index2];

        // Interpolate scalar value at the target point
        double interpolated_value = interpolate_scalar(target, vertex0, vertex1, vertex2, scalar0, scalar1, scalar2);

        // Output the result
        std::cout << "Interpolated scalar at (" << target << ") = " << interpolated_value << "\n";
    }

    return 0;
}