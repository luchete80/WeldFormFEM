#include "Omega_h_r3d.hpp"

static void test_3d() {

  r3d::Few<r3d::Vector<3>, 4> verts = {
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  auto faces = r3d::faces_from_verts(verts);
  OMEGA_H_CHECK(Omega_h::are_close(Omega_h::from_r3d(faces[0].n),
      -Omega_h::normalize(Omega_h::vector_3(1, 1, 1))));
  OMEGA_H_CHECK(Omega_h::are_close(faces[0].d, -faces[0].n[0]));+
  
}
