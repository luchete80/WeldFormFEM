/*************************************************************************/
/*  cylinder_adapt_tst.cpp                                       */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                */
/*  ('https://www.opensourcemech.com',)                                    */
/*                                                                       */
/*  Copyright (c) 2023-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. */ 
/*  See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/




#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto cmdline = Omega_h::CmdLine();
  cmdline.add_arg<std::string>("mesh-path");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto path = cmdline.get<std::string>("mesh-path");
  auto mesh = Omega_h::gmsh::read(path, lib.world());
  auto opts = Omega_h::AdaptOpts(&mesh);
  opts.verbosity = Omega_h::EXTRA_STATS;
  auto metric_input = Omega_h::MetricInput();
  metric_input.should_limit_lengths = true;
  metric_input.max_length = 0.5;
  metric_input.min_length = 0.0;
  metric_input.should_limit_gradation = true;
  metric_input.max_gradation_rate = 1.0;
  auto curvature_source =
      Omega_h::MetricSource(OMEGA_H_CURVATURE, Omega_h::PI / 16.0);
  metric_input.sources.push_back(curvature_source);
  Omega_h::add_implied_metric_tag(&mesh);
  Omega_h::generate_target_metric_tag(&mesh, metric_input);
  while (Omega_h::approach_metric(&mesh, opts)) Omega_h::adapt(&mesh, opts);
  bool ok = check_regression("gold_cylinder", &mesh);
  if (!ok) return 2;
  return 0;
}
