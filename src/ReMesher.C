#include "ReMesher.h"

#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_build.hpp>

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
  ReMesher(Domain_d *d){
    m_dom = d;
  }
  
};