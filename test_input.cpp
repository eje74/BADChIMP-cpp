//#include <vector>
#include "src/io/Input.h"

int main() 
{
    Input input("test.inp");
    input.head_block().info();
    input["version"].info();
    input["dir"].info();
    std::cout << input << std::endl;
    double version = input["version"];                      // version = 1.1
    double pi = input["math"]["pi"];                        // pi = 3.12
    std::vector<int> fibo = input["math"]["series"]["fibo"]; // fibo = 2, 3, 5, 7, 11, 13, 17, 19
    std::vector<int> binary = input["binary"];              // binary = 0, 1, 0, 1, 1, 1, 1, 0
    //int max = input["iterations"]["max"];
    std::cout << version << ", " << pi << std::endl;
    for (const auto& b : binary)
        std::cout << b << ", ";
    std::cout << std::endl;
    std::vector<std::vector<int>> binary2D = input["binary"].matrix<int>();              // binary = 0, 1, 0, 1, 1, 1, 1, 0
    const std::vector<std::string> dir = input["dir"];
    std::cout << dir[0] << std::endl;
}