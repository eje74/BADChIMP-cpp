
//#include <cstdio>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <climits>
#include <cfloat>
#include <string>
#include <iomanip>
#include <map>
#include <sys/types.h>
#if defined (_WIN32)
#include <direct.h>
#else
#include <unistd.h>
#endif
#include <sys/stat.h>
#include <bitset>
#include "Input.h"
#include "Output.h"
#include "MPI.h"

int main(int argc, char *argv[]) {

  Input in;
  in.read("input.dat");  //in.print();

  std::vector<int> n = {10,10,10};

  MPI mpi(&argc, &argv, n, in["mpi"]["procs"]);
  //mpi.print();

  int size = mpi.get_local_total_size();
  std::vector<double> u_tot(n.size()*size);
  std::vector<double> rho(size);
  
  Output output("out", mpi, nullptr);
  output.add_file("fluid");
  output["fluid"].add_variables( {"velocity","density"}, {&u_tot[0],&rho[0]}, {sizeof(u_tot[0]),sizeof(rho[0])}, {3,1}, {3,1});
  output["fluid"].write(0.0);
  //output["fluid"].write(1.0);

  mpi.end();

  //output.add_file("chem");


}


