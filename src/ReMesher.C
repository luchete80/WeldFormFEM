#include "ReMesher.h"
#include "Domain_d.h"

//#include <Omega_h.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_quality.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_timer.hpp>
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
    build_from_elems_and_coords(&mesh,OMEGA_H_SIMPLEX, 3, device_tets, device_coords); // Correct method

    // Step 2: Add the node coordinates as a tag (e.g., "coords")
    mesh.set_tag(Omega_h::VERT, "metric", Omega_h::Reals(device_coords));

    // Step 3: Add element connectivity as a tag (e.g., "conn")
    //mesh.set_tag(Omega_h::CELL, "conn", Omega_h::Reals(device_tets));

    // Optionally, you can also print out to verify the number of nodes and elements
    //std::cout << "Mesh created with " << mesh.nverts() << " vertices and "
    //          << mesh.nelems() << " elements.\n";
}


// Function to evaluate if refinement is necessary based on the scalar field (e.g., plastic strain)
bool needs_refinement(double scalar_value, double threshold) {
    return scalar_value > threshold;
}


/////////////////////////// MANUAL APPROACHES
/*
void refine_mesh(Omega_h::Mesh& mesh, double threshold, 
                 Omega_h::Write<Omega_h::Real>& plastic_strain_field) {
    Omega_h::LOs elem_to_nodes = mesh.ask_down(3, 0).ab2b;  // Get element-node connectivity
    Omega_h::Reals node_coords = mesh.coords();             // Get node positions

    Omega_h::HostWrite<Omega_h::Real> new_coords;           // Storage for new coordinates
    Omega_h::HostWrite<Omega_h::LO> new_conn;               // Storage for new connectivity

    int num_existing_nodes = mesh.nverts();                 // Number of nodes before adding new ones
    int num_existing_elements = mesh.nents(3);              // Number of tetrahedral elements

    for (int elem = 0; elem < num_existing_elements; ++elem) {
        double max_strain = 0.0;
        int nodes[4];

        // Get the nodes of this tetrahedron
        for (int j = 0; j < 4; ++j) {
            nodes[j] = elem_to_nodes[elem * 4 + j];
            double strain = plastic_strain_field[nodes[j]];
            max_strain = std::max(max_strain, strain);
        }

        // If the strain exceeds the threshold, add new nodes
        if (needs_refinement(max_strain, threshold)) {
            std::cout << "Element " << elem << " needs refinement due to strain: " << max_strain << "\n";

            // Compute new node positions (e.g., midpoint of first two nodes)
            int node1 = nodes[0];
            int node2 = nodes[1];

            double x_new = 0.5 * (node_coords[node1 * 3] + node_coords[node2 * 3]);
            double y_new = 0.5 * (node_coords[node1 * 3 + 1] + node_coords[node2 * 3 + 1]);
            double z_new = 0.5 * (node_coords[node1 * 3 + 2] + node_coords[node2 * 3 + 2]);

            // Add the new node to the list
            new_coords.push_back(x_new);
            new_coords.push_back(y_new);
            new_coords.push_back(z_new);

            // Get the new node index
            int new_node_index = num_existing_nodes + new_coords.size() / 3 - 1;

            // Update connectivity by inserting new node into tetrahedron (replacing node 1)
            for (int j = 0; j < 4; ++j) {
                if (j == 1) {
                    new_conn.push_back(new_node_index);
                } else {
                    new_conn.push_back(nodes[j]);
                }
            }
        }
    }

    // Merge new nodes into mesh
    if (new_coords.size() > 0) {
        mesh.add_tag(0, "new_coords", 3, Omega_h::Reals(new_coords));  // Store new node positions
        mesh.add_tag(3, "new_conn", 4, Omega_h::LOs(new_conn));         // Store new connectivity
        std::cout << "Added " << new_coords.size() / 3 << " new nodes.\n";
    }
}
*/
/*
void refine_mesh(Omega_h::Mesh& mesh, double threshold, 
                 Omega_h::Write<Omega_h::Real>& plastic_strain_field) {
    // Create a refinement target
    Omega_h::Reals refine_criteria(mesh.nelems(), 0.0);

    // Mark elements that need refinement
    for (int elem = 0; elem < mesh.nelems(); ++elem) {
        double max_strain = 0.0;

        // Get maximum strain in element
        for (int node = 0; node < 4; ++node) {
            max_strain = std::max(max_strain, plastic_strain_field[mesh.ask_down(3, 0).ab2b[elem * 4 + node]]);
        }

        // Mark for refinement if above threshold
        refine_criteria[elem] = (max_strain > threshold) ? 1.0 : 0.0;
    }

    // Apply refinement operation
    Omega_h::AdaptOpts adapt_opts(&mesh);
    adapt_opts.xfer_opts[Omega_h::OMEGA_H_METRIC].type = Omega_h::OMEGA_H_CONSERVE;
    adapt_opts.refine_fraction = 1.0; // Allow full refinement

    mesh.set_tag(Omega_h::REGION, "refinement", refine_criteria);
    mesh.adapt(adapt_opts);

    std::cout << "Mesh refinement completed based on plastic strain.\n";
}
*/

/*
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
*/
//FROM test_degree

static void tet_run(Omega_h::Real side_angle_in_degrees) {
  auto side_angle = side_angle_in_degrees / 180. * Omega_h::PI;
  std::cout << "side_angle " << side_angle << '\n';
  auto dihedral_angle =
      std::acos((std::cos(side_angle) - Omega_h::square(cos(side_angle))) /
                (Omega_h::square(std::sin(side_angle))));
  std::cout << "dihedral_angle " << dihedral_angle << '\n';
  std::cout << "dihedral_angle in degrees "
            << (dihedral_angle * 180. / Omega_h::PI) << '\n';
  auto solid_angle = 3. * dihedral_angle - Omega_h::PI;
  std::cout << "solid_angle " << solid_angle << '\n';
  auto degree = 4. * Omega_h::PI / solid_angle;
  std::cout << "degree " << degree << '\n';
  auto half_side_angle = side_angle / 2.;
  auto half_cord_length = std::sin(half_side_angle);
  auto cord_length = 2. * half_cord_length;
  std::cout << "cord_length " << cord_length << '\n';
  auto surf_tri_height = std::sqrt(3.) / 2. * cord_length;
  auto inradius = surf_tri_height / 3.;
  auto circumradius = inradius * 2.;
  std::cout << "circumradius " << circumradius << '\n';
  auto surf_tri_area = std::sqrt(3.) / 4. * Omega_h::square(cord_length);
  auto tet_height = std::sqrt(1. - Omega_h::square(circumradius));
  auto tet_volume = 1. / 3. * surf_tri_area * tet_height;
  std::cout << "tet_volume " << tet_volume << '\n';
  auto msl = (3 * Omega_h::square(cord_length) + 3.) / 6.;
  auto quality = Omega_h::mean_ratio<3>(tet_volume, msl);
  std::cout << "quality " << quality << '\n';
}

/*
// Function to refine mesh based on element quality
void refine_mesh_quality(Omega_h::Mesh& mesh, double quality_threshold) {
    // Compute element qualities using a predefined metric (3D element, using the default metric dimension)
    auto quality_measure = Omega_h::MetricElementQualities<3, 3>(&mesh);

    // Create a writable array for refinement criteria
    Omega_h::Write<Omega_h::Real> refine_criteria(mesh.nelems(), 0.0);

    // Identify elements with poor quality
    for (int elem = 0; elem < mesh.nelems(); ++elem) {
        if (quality_measure[elem] < quality_threshold) {
            refine_criteria[elem] = 1.0;  // Mark for refinement
        }
    }

    // Set refinement tag
    if (!mesh.has_tag(Omega_h::REGION, "refinement")) {
        mesh.add_tag<Omega_h::Real>(Omega_h::REGION, "refinement", 1, 0.0);
    }
    mesh.set_tag(Omega_h::REGION, "refinement", refine_criteria);

    // Set up adaptation options
    Omega_h::AdaptOpts adapt_opts(&mesh);
    adapt_opts.verbosity = Omega_h::EXTRA_STATS; // Optional: controls output verbosity

    // Perform h-adaptation (refinement based on element quality)
    Omega_h::adapt(&mesh, adapt_opts);

    std::cout << "Simple quality adaptation completed.\n";
}
*/
using namespace Omega_h;
/*
static void set_target_metric(Mesh* mesh) {
  int dim = 3;
  auto coords = mesh->coords();
  auto target_metrics_w = Write<Real>(mesh->nverts() * symm_ncomps(3));
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto z = coords[v * 3 + (3 - 1)];
    auto h = Vector<3>();
    for (Int i = 0; i < 3 - 1; ++i) h[i] = 0.1;
    h[3 - 1] = 0.001 + 0.198 * std::abs(z - 0.5);
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  parallel_for(mesh->nverts(), f);
  mesh->set_tag(VERT, "target_metric", Reals(target_metrics_w));
}

void refine_mesh_quality(Omega_h::Mesh* mesh, double quality_threshold) {
  int dim = 3;
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = get_implied_metrics(mesh);
  mesh->add_tag(VERT, "metric", symm_ncomps(dim), implied_metrics);
  mesh->add_tag<Real>(VERT, "target_metric", symm_ncomps(dim));
  set_target_metric<dim>(mesh);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->ask_lengths();
  mesh->ask_qualities();
  vtk::FullWriter writer;
  if (vtk_path) {
    writer = vtk::FullWriter(vtk_path, mesh);
    writer.write();
  }
  auto opts = AdaptOpts(mesh);
  opts.verbosity = EXTRA_STATS;
  opts.length_histogram_max = 2.0;
  opts.max_length_allowed = opts.max_length_desired * 2.0;
  Now t0 = now();
  while (approach_metric(mesh, opts)) {
    adapt(mesh, opts);
    if (mesh->has_tag(VERT, "target_metric")) set_target_metric<dim>(mesh);
    if (vtk_path) writer.write();
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";
}
*/

/////// FROM CYLINDER ADAPT_TEST
void refine_mesh_quality(Omega_h::Mesh &mesh){
    // Step 1: Set up adaptation options
    auto opts = Omega_h::AdaptOpts(&mesh);
    opts.verbosity = Omega_h::EXTRA_STATS;  // Optional: set verbosity for more detailed output

    // Step 2: Set up the metric input (mesh quality control)
    auto metric_input = Omega_h::MetricInput();
    metric_input.should_limit_lengths = true;
    metric_input.max_length = 0.5;  // Maximum length for the elements
    metric_input.min_length = 0.0;  // Minimum length for the elements
    metric_input.should_limit_gradation = true;
    metric_input.max_gradation_rate = 1.0;  // Max allowed gradation rate for element sizes

    // Step 3: Set up curvature as a metric source
    //auto curvature_source = Omega_h::MetricSource(OMEGA_H_CURVATURE, Omega_h::PI / 16.0);  // Curvature source with angle tolerance
    //metric_input.sources.push_back(curvature_source);
    //auto solid_angle_source = Omega_h::MetricSource(Omega_h::OMEGA_H_MIN_SOLID_ANGLE);
    //solid_angle_metric_input.sources.push_back(solid_angle_source);


    mesh.set_parting(OMEGA_H_ELEM_BASED);
    mesh.ask_lengths();
    mesh.ask_qualities();

    // Step 4: Generate target metric tag for the mesh
    Omega_h::add_implied_metric_tag(&mesh);  // Add implied metric tag
    Omega_h::generate_target_metric_tag(&mesh, metric_input);  // Generate target metric for mesh

    // Step 5: Approach the metric (ensure mesh adapts to the desired quality)
    while (Omega_h::approach_metric(&mesh, opts)) {
        // Step 6: Adapt the mesh using the metric
        Omega_h::adapt(&mesh, opts);
    }
    
    std::cout << "Mesh refinement based on quality completed.\n";
}


/// TEST THIS TO REPLACE METRIC
/*
   // Define aspect ratio as a metric source
    auto aspect_ratio_source = Omega_h::MetricSource(Omega_h::OMEGA_H_ASPECT_RATIO);
    aspect_ratio_metric_input.sources.push_back(aspect_ratio_source);

    // Add metric tags and generate the target metric tag for elements
    Omega_h::add_implied_metric_tag(&mesh);
    auto target_metrics_w = Omega_h::Write<Omega_h::Real>(mesh.nelems() * 6); // 6 for tetrahedral elements
    Omega_h::generate_target_metric_tag(&mesh, aspect_ratio_metric_input);

    // Set the tag for ELEMENT to store the metric values
    mesh.set_tag(Omega_h::ELEMENT, "target_metric", target_metrics_w);

    // Perform mesh adaptation based on the chosen metric
    while (Omega_h::approach_metric(&mesh, opts)) {
        Omega_h::adapt(&mesh, opts);
    }
*/

namespace MetFEM{
  ReMesher::ReMesher(Domain_d *d){
    m_dom = d;


    Omega_h::Library lib;
    Omega_h::Mesh mesh(&lib);

    // Example data
    /*
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
    
    */

#ifdef CUDA_BUILD
    // Allocate GPU memory
    double* d_node_data;
    int* d_connectivity_data;
    cudaMalloc((void**)&d_node_data, m_dom->m_node_count * 3 * sizeof(double));
    cudaMalloc((void**)&d_connectivity_data, num_elements * 4 * sizeof(int));

    // Copy from host to GPU
    cudaMemcpy(d_node_data, x, num_nodes * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_connectivity_data, h_connectivity_data, num_elements * 4 * sizeof(int), cudaMemcpyHostToDevice);

    // Create mesh using GPU data
    create_mesh(mesh, d_node_data, num_nodes, d_connectivity_data, num_elements);

    // Free GPU memory
    cudaFree(d_node_data);
    cudaFree(d_connectivity_data);
#else
    // CPU case
    create_mesh(mesh, m_dom->x, m_dom->m_node_count, (int *)m_dom->m_elnod, m_dom->m_elem_count);
#endif
    std::cout << "Generating remeshing"<<std::endl;
    refine_mesh_quality(mesh);
    std::cout << "Refine done "<<std::endl;
    // Save mesh
    //Omega_h::write_mesh("output.osh", &mesh);
    //Omega_h::vtk_export_mesh("output.vtk", &mesh);
    Omega_h::vtk::Writer writer("out_amr_3D", &mesh);
    auto w = lib.world();
    writer.write();
    
    

  // auto lib_osh = Omega_h::Library(&argc, &argv);
  // auto comm_osh = lib_osh.world();

  // auto mesh_osh = Omega_h::build_box(comm_osh,
      // OMEGA_H_SIMPLEX,
      // 1.0, 1.0, 0.0, 32, 32, 0);
  // mesh_osh.balance();
  // mesh_osh.set_parting(OMEGA_H_GHOSTED);
  // Omega_h::add_implied_metric_tag(&mesh_osh);
  // mesh_osh.set_parting(OMEGA_H_ELEM_BASED);

// #ifdef OMEGA_H_USE_MPI
  // auto mesh_dolfin = std::make_shared<dolfin::Mesh>(comm_osh->get_impl());
// #else
  // auto mesh_dolfin = std::make_shared<dolfin::Mesh>();
// #endif

  // dolfin::File file_dolfin("dolfin.pvd");
  // Omega_h::vtk::Writer file_osh("omega_h-vtk-output", &mesh_osh);

  // int i = 0;
  // int n = 3;
  // while (true) {

    // Omega_h::to_dolfin(*mesh_dolfin, &mesh_osh);

    // auto V = std::make_shared<Poisson::FunctionSpace>(mesh_dolfin);

    // auto u0 = std::make_shared<dolfin::Constant>(0.0);
    // auto boundary = std::make_shared<DirichletBoundary>();
    // dolfin::DirichletBC bc(V, u0, boundary);

    // Poisson::BilinearForm a(V, V);
    // Poisson::LinearForm L(V);
    // auto f = std::make_shared<Source>();
    // auto g = std::make_shared<dUdN>();
    // L.f = f;
    // L.g = g;

    // dolfin::Function u(V);
    // solve(a == L, u, bc);

    // Omega_h::from_dolfin(&mesh_osh, u, "u");

    // file_dolfin << u;
    // file_osh.write();

    // if (++i == n) break;

    // mesh_osh.set_parting(OMEGA_H_GHOSTED);
    // Omega_h::MetricInput metric_input;
    // auto source = Omega_h::MetricSource(OMEGA_H_VARIATION, 2e-3, "u");
    // metric_input.sources.push_back(source);
    // metric_input.should_limit_lengths = true;
    // metric_input.max_length = 1.0 / 2.0;
    // metric_input.should_limit_gradation = true;
    // Omega_h::generate_target_metric_tag(&mesh_osh, metric_input);
    // Omega_h::AdaptOpts opts(&mesh_osh);
    // opts.verbosity = Omega_h::EXTRA_STATS;
    // while (Omega_h::approach_metric(&mesh_osh, opts)) {
      // Omega_h::adapt(&mesh_osh, opts);
    // }

  // }


  }
  
};
