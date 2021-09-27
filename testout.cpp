#include <vector>
#include <array>
#include <iostream>
#include "src/io/vector_func.h"
#include "src/io/VTK.h"

int main()
{
    // std::vector<std::vector<double>> nodes = {{0,0,0},{0,1,0},{0,0,1},{0,1,1}};
    // std::vector<int> size = {2,2,2};
    //std::vector<std::vector<float>> nodes = {{1,1,1},{4,4,4}};
    std::vector<std::vector<float>> nodes = {{1,1},{1,2},{1,3}, {2,1},{2,2},{2,3}, {3,1},{3,2},{3,3}};
    VTK::Output<VTK::pixel, double> out(nodes, "out_test", 0, 1, VTK::ASCII); 
    out.add_file("fluid");
    auto data = linspace(1.0, 10.0, 1.0);
    //std::cout << data << std::endl;
    out.add_variable("rho", 1, data);
    out.write(0);
}
