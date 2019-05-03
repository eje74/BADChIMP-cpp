#include "Input.h"
#include "Output.h"
#include "MPI.h"


int main(int argc, char *argv[]) {

  // read input files
  Input param("input.dat"); //param.print();
  Input geo("SK_geo.dat"); //geo.print();

  // initialize MPI
  MPI mpi(&argc, &argv, geo["dim"], param["mpi"]["procs"]);
  mpi.print();

  // allocate and initialize geometry
  int size = mpi.get_local_size();
  std::vector<char> nodes(size);
  mpi.set_geometry(nodes, geo["geo"]);

  // allocate macroscopic variables
  std::vector<double> u_tot(mpi.get_dim()*size);
  std::vector<double> rho(size);
  
  // initialize output
  Output output("out", mpi, nullptr);
  output.add_file("geo");
  output["geo"].add_variables({"geo"}, {&nodes[0]}, {sizeof(nodes[0])}, {1}, {1});
  output.add_file("fluid");
  output["fluid"].add_variables({"velocity","density"}, {&u_tot[0],&rho[0]}, {sizeof(u_tot[0]),sizeof(rho[0])}, {3,1}, {3,1});

  // write output
  output["geo"].write(0);
  output["fluid"].write(0);

  mpi.end();

}
