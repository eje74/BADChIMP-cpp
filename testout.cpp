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
    auto rho_data = linspace(1.0, 20.0, 1.0);
    //std::valarray<double> arr_data(data.data(), data.size());
    auto index1 = linspace(1, 18+3, 3);
    auto index2 = linspace(0, 16+3, 3);
    //std::vector<double> vel = { 0.1,0.1,0.1, 0.2,0.2,0.2, 0.3,0.3,0.3, 0.4,0.4,0.4, 0.5,0.5,0.5, 0.6,0.6,0.6, 0.7,0.7,0.7, 0.8,0.8,0.8, 0.9,0.9,0.9 };
    // auto rho_buf = out.add_buffer(arr_data);
    // out.add_variable("rho1", 1, rho_buf.data(), index1);
    // out.add_variable("rho2", 1, rho_buf.data(), index2);
    out.add_variable("rho1", 1, rho_data, index1);
    out.add_variable("rho2", 1, rho_data, index2); 
   // out.add_variable("rho1", 1, vec_data, index1);
    // out.add_variable("rho2", 1, data, index2);
    std::valarray<double> vel = { 0.1,0.1,0.1, 0.2,0.2,0.2, 0.3,0.3,0.3, 0.4,0.4,0.4, 0.5,0.5,0.5, 0.6,0.6,0.6, 0.7,0.7,0.7, 0.8,0.8,0.8, 0.9,0.9,0.9 };
    //auto& vel_buf = out.add_buffer(vel);
    out.add_buffer(vel);
    //std::cout << " test: " << vel_buf.data().data() << " " << vel_buf.data() << std::endl;
    //vel_buf.update();
    // out.add_variable("vel", 3, vel_buf.data());
    out.add_variable("vel", 3, out.buffers().back().data());
    out.write(0);
}
