#include "ReMesher.h"
#include "Domain_d.h"

//#include <Omega_h.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_amr.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_quality.hpp>
#include <Omega_h_map.hpp> //colloect matrked
#include <Omega_h_metric.hpp>
#include <Omega_h_timer.hpp>
#include <Omega_h_coarsen.hpp>
#include "Omega_h_refine_qualities.hpp"
#include "Omega_h_array_ops.hpp"      /// ARE CLOSE
#include "Omega_h_class.hpp"

#ifdef CUDA_BUILD
#include <cuda_runtime.h>
#endif

void classify_vertices(Omega_h::Mesh* mesh) {
    // Assume all vertices are in the interior (volume = 3)
    Omega_h::Write<Omega_h::I8> class_dim(mesh->nverts(), 3);

    // If you have known boundary conditions, modify class_dim accordingly
    // Example: If you have a function that detects surface nodes, set them to 2

    // Attach classification data to mesh
    mesh->add_tag(Omega_h::VERT, "class_dim", 1, Omega_h::Read<Omega_h::I8>(class_dim));
}

/*
void classify_vertices(Omega_h::Mesh* mesh) {
    // Assume all vertices are in the interior (volume = 3)
    Omega_h::Write<Omega_h::I8> class_dim(mesh->nverts(), 3); 

    // Example: Classify boundary vertices
    auto coords = mesh->coords();
    auto f = OMEGA_H_LAMBDA(Omega_h::LO v) {
        auto x = get_vector<3>(coords, v);
        // Check if the vertex is near the boundary (e.g., within a small tolerance)
        if (x[0] < 0.01 || x[1] < 0.01 || x[2] < 0.01) {
            class_dim[v] = 2; // Boundary classification
        } else {
            class_dim[v] = 3; // Interior classification
        }
    };
    parallel_for(mesh->nverts(), f);

    // Attach classification to the mesh
    mesh->add_tag(Omega_h::VERT, "class_dim", 1, Omega_h::Read<Omega_h::I8>(class_dim));
}
*/

using namespace Omega_h;
template <Int dim>
void create_class_dim(Mesh* mesh) {
    // Number of vertices in the mesh
    auto nverts = mesh->nverts();
    
    // Retrieve the coordinates of the vertices
    auto coords = mesh->coords();
    
    // Create a Write<Byte> container for storing the classification results
    auto class_dims = Write<Byte>(nverts);

    // Define some boundary criteria (for example, 0.5 could represent a boundary threshold)
    Real boundary_threshold = 0.5;

    // Loop over all the vertices to classify them
    auto f = OMEGA_H_LAMBDA(LO v) {
        // Get the coordinates of the current vertex (assuming 3D mesh here)
        auto x = get_vector<dim>(coords, v);
        
        // Simple classification logic based on some threshold (e.g., x[0] > boundary_threshold)
        Byte class_dim = 0; // Initially set to "interior"
        
        // For example, if the vertex is close to the boundary, we classify it as boundary
        if (x[0] > boundary_threshold) {
            class_dim = 1;  // "boundary"
        }
        
        // Store the classification result in the Write<Byte> container
        class_dims[v] = class_dim;
    };

    // Run the classification function in parallel
    parallel_for(nverts, f, "classify_vertices");

    // Add the class_dim tag to the mesh
    mesh->add_tag<Byte>(VERT, "class_dim", 1, class_dims);
}


template <Int dim>
void classify_faces(Mesh* mesh) {
    // Number of faces in the mesh
    auto nfaces = mesh->nfaces();

    // Retrieve the coordinates of the faces or the centroids of faces
    //auto centroids = mesh->centroids();

    // Create a Write<Byte> container to store the classification results
    auto class_dims = Write<Byte>(nfaces,0);

    // Define some classification logic, for example, boundary faces are classified as 1
    Real boundary_threshold = 0.5;
    /*
    // Loop over all the faces (triangles, quadrilaterals, etc.)
    auto f = OMEGA_H_LAMBDA(LO f) {
        // Get the centroid of the current face (assuming 3D mesh here)
        auto x = get_vector<dim>(centroids, f);
        
        // Initialize classification to interior (0)
        Byte class_dim = 0;

        // Example condition: If the centroid is close to the boundary, classify it as boundary (1)
        if (x[0] > boundary_threshold) {
            class_dim = 1; // Boundary
        }

        // Store the classification result in the Write<Byte> container
        class_dims[f] = class_dim;
    };
    */

    // Run the classification function in parallel
    //parallel_for(nfaces, f, "classify_faces");

    // Add the class_dim tag to the mesh for faces
    mesh->add_tag<Byte>(FACE, "class_dim", 1, class_dims);
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
    std::cout << "Creating nodes "<<std::endl;
    // CPU Case: Copy raw pointer data to Omega_h::HostWrite<>
    Omega_h::HostWrite<Omega_h::Real> coords(num_nodes * 3);
    Omega_h::HostWrite<Omega_h::LO> tets(num_elements * 4);

    std::cout << "Done "<<std::endl;    
    for (int i = 0; i < num_nodes * 3; ++i) coords[i] = h_node_coords[i];
    for (int i = 0; i < num_elements * 4; ++i) tets[i] = h_element_conn[i];

    std::cout << "Convert to write "<<std::endl;    
    // Convert HostWrite to Write<> for Omega_h
    Omega_h::Write<Omega_h::Real> device_coords(coords);
    Omega_h::Write<Omega_h::LO> device_tets(tets);
#endif
    std::cout << "Building from elements "<<std::endl;
    // Build mesh (works on both CPU and GPU)
    build_from_elems_and_coords(&mesh,OMEGA_H_SIMPLEX, 3, device_tets, device_coords); // Correct method

  if (!mesh.has_tag(Omega_h::VERT, "coordinates")) {
      std::cerr << "Error: Mesh does not have 'coordinates' tag!" << std::endl;
  }
  classify_elements(&mesh);
  classify_faces<3>(&mesh);
  classify_vertices(&mesh);
  create_class_dim<3>(&mesh);
  
if (!mesh.has_tag(Omega_h::VERT, "class_dim")) {
    std::cerr << "Error: Mesh does not have 'class_dim' tag!" << std::endl;
}
    // Step 2: Add the node coordinates as a tag (e.g., "coords")
    //mesh.set_tag(Omega_h::VERT, "metric", Omega_h::Reals(device_coords));

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
//// GIVES assertion false failed at /home/weldform-pc/Numerico/WeldFormFEM/lib/omega_h-9.34.13/src/Omega_h_refine_qualities.cpp +106

/////// FROM UNIT_MESH (several tests)
void refine_mesh_quality(Omega_h::Mesh &mesh){
  //build_box_internal(&mesh, OMEGA_H_SIMPLEX, 1., 1., 0., 1, 1, 0);
  LOs candidates(mesh.nedges(), 0, 1);
  mesh.add_tag(VERT, "metric", symm_ncomps(2),
      repeat_symm(mesh.nverts(), identity_matrix<2, 2>()));
  auto quals = refine_qualities(&mesh, candidates);
  OMEGA_H_CHECK(are_close(
      quals, Reals({0.494872, 0.494872, 0.866025, 0.494872, 0.494872}), 1e-4));
}
/*
void adapt (Omega_h::Mesh &mesh) {
  
  auto opts = AdaptOpts(&mesh);
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(0.3)));
  while (coarsen_by_size(&mesh, opts))
    ;
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(0.6)));
  while (coarsen_by_size(&mesh, opts))
    ;
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(1.0)));
  while (coarsen_by_size(&mesh, opts))
    ;
  mesh.ask_qualities();
  bool ok = check_regression("gold_coarsen", &mesh);  
  
}
*/
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


//////////////////////////// FROM UGAWG LINEAR


template <Int dim>
static void set_target_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto target_metrics_w = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto z = coords[v * dim + (dim - 1)];
    auto h = Vector<dim>();
    for (Int i = 0; i < dim - 1; ++i) h[i] = 0.1;
    h[dim - 1] = 0.001 + 0.198 * std::abs(z - 0.5);
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  parallel_for(mesh->nverts(), f);
  mesh->set_tag(VERT, "target_metric", Reals(target_metrics_w));
}

template <Int dim>
void run_case(Mesh* mesh, char const* vtk_path) {
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = get_implied_metrics(mesh);
  mesh->add_tag(VERT, "metric", symm_ncomps(dim), implied_metrics);
  std::cout << "symm_ncomps(" << dim << ") = " << symm_ncomps(dim) << std::endl;
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
  std::cout << "Adapting "<<std::endl; 
  while (approach_metric(mesh, opts)) {
    std::cout << "Step "<<std::endl;
    adapt(mesh, opts);
    std::cout << "DONE"<<std::endl;
    //if (mesh->has_tag(VERT, "target_metric")) set_target_metric<dim>(mesh);
    //else      std::cerr << "Error: target_metric tag was not properly set!" << std::endl;
    //if (vtk_path) writer.write();
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";
}

////////////////////////////////////////////////////////////////////////////////////////////
/// FROM AMR_TEST2

template <int dim>
OMEGA_H_INLINE double eval_rc(Omega_h::Vector<dim> c);

template <>
OMEGA_H_INLINE double eval_rc<2>(Omega_h::Vector<2> c) {
  auto rc2 = (c[0] - 0.5) * (c[0] - 0.5) + (c[1] - 0.5) * (c[1] - 0.5);
  return std::sqrt(rc2);
}

template <>
OMEGA_H_INLINE double eval_rc<3>(Omega_h::Vector<3> c) {
  auto rc2 = (c[0] - 0.5) * (c[0] - 0.5) + (c[1] - 0.5) * (c[1] - 0.5) +
             (c[2] - 0.5) * (c[2] - 0.5);
  return std::sqrt(rc2);
}


template <int dim>
Omega_h::Bytes mark(Omega_h::Mesh* m, int level) {
  auto coords = m->coords();
  auto mids = Omega_h::average_field(m, dim, dim, coords);
  auto is_leaf = m->ask_leaves(dim);
  auto leaf_elems = Omega_h::collect_marked(is_leaf);
  Omega_h::Write<Omega_h::Byte> marks(m->nelems(), 0);
  auto f = OMEGA_H_LAMBDA(Omega_h::LO e) {
    auto elem = leaf_elems[e];
    auto c = Omega_h::get_vector<dim, Omega_h::Reals>(mids, elem);
    auto rc = eval_rc<dim>(c);
    auto r = 0.25;
    auto tol = 0.314 / static_cast<double>(level);
    if (std::abs(rc - r) < tol) marks[elem] = 1;
  };
  Omega_h::parallel_for(leaf_elems.size(), f);
  return Omega_h::amr::enforce_2to1_refine(m, dim - 1, marks);
}


void refine(Mesh* mesh){
  Omega_h::vtk::Writer writer("out_amr_3D", mesh);
  writer.write();
  for (int i = 1; i < 5; ++i) {
    auto xfer_opts = Omega_h::TransferOpts();
    auto marks = mark<3>(mesh, i);
    Omega_h::amr::refine(mesh, marks, xfer_opts);
    writer.write();
  }
}


//////////////////////////////////////////////////////////////////////////////////////
// FOM WARP

using namespace Omega_h;

void add_dye(Mesh* mesh) {
  auto dye_w = Write<Real>(mesh->nverts());
  auto coords = mesh->coords();
  auto dye_fun = OMEGA_H_LAMBDA(LO vert) {
    auto x = get_vector<3>(coords, vert);
    auto left_cen = vector_3(.25, .5, .5);
    auto right_cen = vector_3(.75, .5, .5);
    auto left_dist = norm(x - left_cen);
    auto right_dist = norm(x - right_cen);
    auto dist = min2(left_dist, right_dist);
    if (dist < .25) {
      auto dir = sign(left_dist - right_dist);
      dye_w[vert] = 4.0 * dir * (.25 - dist);
    } else {
      dye_w[vert] = 0;
    }
  };
  parallel_for(mesh->nverts(), dye_fun);
  mesh->add_tag(VERT, "dye", 1, Reals(dye_w));
}

Reals form_pointwise(Mesh* mesh) {
  auto dim = mesh->dim();
  auto ecoords =
      average_field(mesh, dim, LOs(mesh->nelems(), 0, 1), dim, mesh->coords());
  auto pw_w = Write<Real>(mesh->nelems());
  auto pw_fun = OMEGA_H_LAMBDA(LO elem) { pw_w[elem] = ecoords[elem * dim]; };
  parallel_for(mesh->nelems(), pw_fun);
  return pw_w;
}

static void add_pointwise(Mesh* mesh) {
  auto data = form_pointwise(mesh);
  mesh->add_tag(mesh->dim(), "pointwise", 1, data);
}

static void check_total_mass(Mesh* mesh) {
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  Reals masses = multiply_each(densities, sizes);
  auto owned_masses = mesh->owned_array(mesh->dim(), masses, 1);
  OMEGA_H_CHECK(are_close(1.0, get_sum(mesh->comm(), owned_masses)));
}

static void postprocess_pointwise(Mesh* mesh) {
  auto data = mesh->get_array<Real>(mesh->dim(), "pointwise");
  auto expected = form_pointwise(mesh);
  auto diff = subtract_each(data, expected);
  mesh->add_tag(mesh->dim(), "pointwise_err", 1, diff);
}

template <int dim>
void adapt_warp(Mesh &mesh){

 mesh.set_parting(OMEGA_H_GHOSTED);
  auto metrics = get_implied_isos(&mesh);
  mesh.add_tag(VERT, "metric", 1, metrics);
  add_dye(&mesh);
  mesh.add_tag(mesh.dim(), "density", 1, Reals(mesh.nelems(), 1.0));
  add_pointwise(&mesh);
  auto opts = AdaptOpts(&mesh);
  opts.xfer_opts.type_map["density"] = OMEGA_H_CONSERVE;
  opts.xfer_opts.integral_map["density"] = "mass";
  opts.xfer_opts.type_map["pointwise"] = OMEGA_H_POINTWISE;
  opts.xfer_opts.type_map["dye"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.integral_diffuse_map["mass"] = VarCompareOpts::none();
  opts.verbosity = EXTRA_STATS;
  auto mid = zero_vector<dim>();
  mid[0] = mid[1] = .5;
  Now t0 = now();
  for (Int i = 0; i < 8; ++i) {
    auto coords = mesh.coords();
    Write<Real> warp_w(mesh.nverts() * dim);
    auto warp_fun = OMEGA_H_LAMBDA(LO vert) {
      auto x0 = get_vector<3>(coords, vert);
      auto x1 = zero_vector<3>();
      x1[0] = x0[0];
      x1[1] = x0[1];
      auto x2 = x1 - mid;
      auto polar_a = std::atan2(x2[1], x2[0]);
      auto polar_r = norm(x2);
      Real rot_a = 0;
      if (polar_r < 0.5) {
        rot_a = (PI / 8) * (2.0 * (0.5 - polar_r));
        if (i >= 4) rot_a = -rot_a;
      }
      auto dest_a = polar_a + rot_a;
      auto dst = x0;
      dst[0] = std::cos(dest_a) * polar_r;
      dst[1] = std::sin(dest_a) * polar_r;
      dst = dst + mid;
      auto w = dst - x0;
      set_vector<3>(warp_w, vert, w);
    };
    parallel_for(mesh.nverts(), warp_fun);
    mesh.add_tag(VERT, "warp", dim, Reals(warp_w));
    while (warp_to_limit(&mesh, opts)) {
      adapt(&mesh, opts);
    }
  }
  Now t1 = now();
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  if (mesh.comm()->rank() == 0) {
    std::cout << "test took " << (t1 - t0) << " seconds\n";
  }
  check_total_mass(&mesh);
  postprocess_pointwise(&mesh);
  bool ok = check_regression("gold_warp", &mesh);  
}

template <int dim>
void adapt_warp_with_threshold(Mesh &mesh, Real length_threshold, Real angle_threshold) {
    mesh.set_parting(OMEGA_H_GHOSTED);
    auto metrics = get_implied_isos(&mesh);
    mesh.add_tag(VERT, "metric", 1, metrics);
    add_dye(&mesh);
    mesh.add_tag(mesh.dim(), "density", 1, Reals(mesh.nelems(), 1.0));
    add_pointwise(&mesh);
    auto opts = AdaptOpts(&mesh);
    opts.xfer_opts.type_map["density"] = OMEGA_H_CONSERVE;
    opts.xfer_opts.integral_map["density"] = "mass";
    opts.xfer_opts.type_map["pointwise"] = OMEGA_H_POINTWISE;
    opts.xfer_opts.type_map["dye"] = OMEGA_H_LINEAR_INTERP;
    opts.xfer_opts.integral_diffuse_map["mass"] = VarCompareOpts::none();
    opts.verbosity = EXTRA_STATS;

    auto mid = zero_vector<dim>();
    mid[0] = mid[1] = .5;
    Now t0 = now();

    for (Int i = 0; i < 8; ++i) {
        auto coords = mesh.coords();
        Write<Real> warp_w(mesh.nverts() * dim);
        auto warp_fun = OMEGA_H_LAMBDA(LO vert) {
            auto x0 = get_vector<3>(coords, vert);
            auto x1 = zero_vector<3>();
            x1[0] = x0[0];
            x1[1] = x0[1];
            auto x2 = x1 - mid;
            auto polar_a = std::atan2(x2[1], x2[0]);
            auto polar_r = norm(x2);
            Real rot_a = 0;
            if (polar_r < 0.5) {
                rot_a = (PI / 8) * (2.0 * (0.5 - polar_r));
                if (i >= 4) rot_a = -rot_a;
            }
            auto dest_a = polar_a + rot_a;
            auto dst = x0;
            dst[0] = std::cos(dest_a) * polar_r;
            dst[1] = std::sin(dest_a) * polar_r;
            dst = dst + mid;
            auto w = dst - x0;
            set_vector<3>(warp_w, vert, w);
        };
        parallel_for(mesh.nverts(), warp_fun);
        mesh.add_tag(VERT, "warp", dim, Reals(warp_w));

        // Compute element edge lengths
        auto elems2verts = mesh.ask_down(dim, VERT);
        auto vert_coords = mesh.coords();
        bool refine_needed = false;
        
        for (LO elem = 0; elem < mesh.nelems(); ++elem) {
            for (int j = 0; j < dim + 1; ++j) {
                for (int k = j + 1; k < dim + 1; ++k) {
                    auto vj = elems2verts.ab2b[elem * (dim + 1) + j];
                    auto vk = elems2verts.ab2b[elem * (dim + 1) + k];
                    auto pj = get_vector<dim>(vert_coords, vj);
                    auto pk = get_vector<dim>(vert_coords, vk);
                    auto edge_length = norm(pk - pj);
                    if (edge_length > length_threshold) {
                        refine_needed = true;
                        break;
                    }
                }
                if (refine_needed) break;
            }
            if (refine_needed) break;
        }
        
        if (refine_needed) {
            while (warp_to_limit(&mesh, opts)) {
                adapt(&mesh, opts);
            }
        }
    }
    Now t1 = now();
    mesh.set_parting(OMEGA_H_ELEM_BASED);
    if (mesh.comm()->rank() == 0) {
        std::cout << "test took " << (t1 - t0) << " seconds\n";
    }
    check_total_mass(&mesh);
    postprocess_pointwise(&mesh);
    bool ok = check_regression("gold_warp", &mesh);
}


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
    
    //vtk::Reader vtk_reader("piece_0.vtu");
    //classify_vertices(&mesh);
    //vtk_reader.read(mesh);

    // Free GPU memory
    cudaFree(d_node_data);
    cudaFree(d_connectivity_data);
#else
  std::cout << "Creating mesh"<<std::endl;
    // CPU case
    create_mesh(mesh, m_dom->x, m_dom->m_node_count, (int *)m_dom->m_elnod, m_dom->m_elem_count);
#endif
    
  
    //refine_mesh_quality(mesh);
    std::cout << "Refine done "<<std::endl;
    // Save mesh
    //Omega_h::write_mesh("output.osh", &mesh);
    //Omega_h::vtk_export_mesh("output.vtk", &mesh);
    Omega_h::vtk::Writer writer("out_amr_3D", &mesh);
    auto w = lib.world();
    writer.write();
    
    double length_tres = 0.005;
    double ang_tres = 0.1;
    //run_case<3>(&mesh, "test");
    //refine(&mesh);
    adapt_warp<3>(mesh);
    writer.write();    
    
    adapt_warp_with_threshold<3>(mesh,length_tres, ang_tres);
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
