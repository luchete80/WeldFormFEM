#include "ReMesher.h"

#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_build.hpp>

#include "Domain_d.h"

//From build

// void build_from_elems_and_coords(
    // Mesh* mesh, Omega_h_Family family, Int edim, LOs ev2v, Reals coords) {
  // auto nverts = coords.size() / edim;
  // build_from_elems2verts(mesh, family, edim, ev2v, nverts);
  // mesh->add_coords(coords);
// }

// void build_box_internal(Mesh* mesh, Omega_h_Family family, Real x, Real y,
    // Real z, LO nx, LO ny, LO nz, bool symmetric) {
  // OMEGA_H_CHECK(nx > 0);
  // OMEGA_H_CHECK(ny >= 0);
  // OMEGA_H_CHECK(nz >= 0);
  // if (ny == 0) {
    // LOs ev2v;
    // Reals coords;
    // make_1d_box(x, nx, &ev2v, &coords);
    // build_from_elems_and_coords(mesh, family, EDGE, ev2v, coords);
  // } else if (nz == 0) {
    // LOs fv2v;
    // Reals coords;
    // make_2d_box(x, y, nx, ny, &fv2v, &coords);
    // if (family == OMEGA_H_SIMPLEX && (!symmetric)) fv2v = tris_from_quads(fv2v);
    // auto fam2 = symmetric ? OMEGA_H_HYPERCUBE : family;
    // build_from_elems_and_coords(mesh, fam2, FACE, fv2v, coords);
    // if (family == OMEGA_H_SIMPLEX && symmetric) tris_from_quads_symmetric(mesh);
  // } else {
    // LOs rv2v;
    // Reals coords;
    // make_3d_box(x, y, z, nx, ny, nz, &rv2v, &coords);
    // if (family == OMEGA_H_SIMPLEX && (!symmetric)) rv2v = tets_from_hexes(rv2v);
    // auto fam2 = symmetric ? OMEGA_H_HYPERCUBE : family;
    // build_from_elems_and_coords(mesh, fam2, REGION, rv2v, coords);
    // if (family == OMEGA_H_SIMPLEX && symmetric) tets_from_hexes_symmetric(mesh);
  // }
// }

// Mesh build_box(CommPtr comm, Omega_h_Family family, Real x, Real y, Real z,
    // LO nx, LO ny, LO nz, bool symmetric) {
  // auto lib = comm->library();
  // auto mesh = Mesh(lib);
  // if (comm->rank() == 0) {
    // build_box_internal(&mesh, family, x, y, z, nx, ny, nz, symmetric);
    // reorder_by_hilbert(&mesh);
    // classify_box(&mesh, x, y, z, nx, ny, nz);
    // mesh.class_sets = get_box_class_sets(mesh.dim());
  // }
  // mesh.set_comm(comm);
  // mesh.balance();
  // return mesh;
// }

namespace MetFEM{
  ReMesher::ReMesher(Domain_d *d){
    m_dom = d;



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