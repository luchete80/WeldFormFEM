#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("input.osh");
  cmdline.add_arg<std::string>("output.exo");
  cmdline.add_flag("-v", "verbose");
  auto& cwflag = cmdline.add_flag(
      "--classify-with", "comma separated \"node_sets,side_sets\"");
  cwflag.add_arg<std::string>("set-types");
  if (!cmdline.parse(comm, &argc, argv) ||
      !Omega_h::CmdLine::check_empty(comm, argc, argv)) {
    cmdline.show_help(comm, argv);
    return -1;
  }
  auto inpath = cmdline.get<std::string>("input.osh");
  auto outpath = cmdline.get<std::string>("output.exo");
  auto verbose = cmdline.parsed("-v");
  int classify_with;
  if (cmdline.parsed("--classify-with")) {
    auto set_types = cmdline.get<std::string>("--classify-with", "set-types");
    classify_with = 0;
    if (set_types.find("node_sets") != std::string::npos) {
      classify_with |= Omega_h::exodus::NODE_SETS;
    }
    if (set_types.find("side_sets") != std::string::npos) {
      classify_with |= Omega_h::exodus::SIDE_SETS;
    }
  } else {
    classify_with = Omega_h::exodus::NODE_SETS | Omega_h::exodus::SIDE_SETS;
  }
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(inpath, lib.world(), &mesh);
  Omega_h::exodus::write(outpath, &mesh, verbose, classify_with);
  return 0;
}
