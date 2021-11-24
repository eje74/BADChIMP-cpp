#include <vector>
#include <valarray>
#include <iostream>
#include "src/io/vector_func.h"
#include "src/io/VTK.h"


// class Out : public VTK::Output<VTK::voxel, double>
// {
//     public:
//     Out(std::vector<std::vector<float>>& nodes, const int rank, const int num_procs)
//         : VTK::Output<VTK::voxel, double>(VTK::BINARY, nodes, "out_test", rank, num_procs) { }


// };

int main()
{
    std::vector<std::vector<float>> nodes = {{1,1,1},{1,2,1},{1,3,1}, {2,1,2},{2,2,2},{2,3,2} };
    VTK::Output<VTK::voxel, double> out(VTK::BINARY, nodes, "out_test", 0, 1); 
    out.add_file("fluid");

    auto data = linspace(1.0, 20.0, 1.0);
    std::valarray<double> rho_data(data.data(), data.size());
    //auto rho_data = linspace(1.0, 20.0, 1.0);

    auto index1 = linspace(1, 18+3, 3);
    auto index2 = linspace(0, 16+3, 3);
    out.add_variable("rho1", 1, rho_data, index1);
    out.add_variable("rho2", 1, rho_data, index2); 

    std::valarray<double> vel = { 0.1,0.1,0.1, 0.2,0.2,0.2, 0.3,0.3,0.3, 0.4,0.4,0.4, 0.5,0.5,0.5, 0.6,0.6,0.6 };
    //std::valarray<double> vel2 = { 1.1,1.1,1.1, 1.2,1.2,1.2, 1.3,1.3,1.3, 1.4,1.4,1.4, 1.5,1.5,1.5, 1.6,1.6,1.6 };
    std::vector<double> vel2 = { 1.1,1.1,1.1, 1.2,1.2,1.2, 1.3,1.3,1.3, 1.4,1.4,1.4, 1.5,1.5,1.5, 1.6,1.6,1.6 };
    out.add_variable("vel", 3, vel);
    out.add_variable("vel2", 3, vel2, linspace(0, (int)vel2.size(), 1));
    out.write(0);
}
