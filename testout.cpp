#include <vector>
#include <array>
#include <iostream>
#include "src/io/vector_func.h"
#include "src/io/VTK.h"

int main()
{
    //std::vector<std::vector<float>> nodes = {{1,1,1},{4,4,4}};
    //std::vector<std::vector<float>> nodes = {{1,1},{1,2},{1,3}, {2,1},{2,2},{2,3}, {3,1},{3,2},{3,3}};
    std::vector<std::vector<float>> nodes = {{1,1,1},{1,2,1},{1,3,1}, {2,1,2},{2,2,2},{2,3,2} };
    //std::vector<std::vector<float>> nodes = {{1,1},{1,2},{1,3}, {2,1},{2,2},{2,3}, {3,1},{3,2},{3,3}};
    //std::vector<std::vector<float>> nodes = {{1,1,1},{1,2,1},{1,3,1}, {2,1,1},{2,2,1},{2,3,1}, {3,1,1},{3,2,1},{3,3,1}};
    VTK::Output<VTK::voxel, double> out(nodes, "out_test", 0, 1, VTK::BINARY); 
    out.add_file("fluid");
    auto data = linspace(1.0, 20.0, 1.0);
    auto index1 = linspace(1, 18, 3);
    auto index2 = linspace(0, 16, 3);
    //std::vector<double> vel = { 0.1,0.1,0.1, 0.2,0.2,0.2, 0.3,0.3,0.3, 0.4,0.4,0.4, 0.5,0.5,0.5, 0.6,0.6,0.6, 0.7,0.7,0.7, 0.8,0.8,0.8, 0.9,0.9,0.9 };
    //std::cout << index << std::endl;
    out.add_variable("rho1", 1, data, index1);
    out.add_variable("rho2", 1, data, index2);
    //out.add_variable("vel", 3, vel);
    //out.add_variable("rho", 1, data);
    out.write(0);
}
