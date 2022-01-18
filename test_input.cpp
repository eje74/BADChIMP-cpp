//#include <vector>
#include "src/io/Input.h"

int main() 
{
    Input input("test.inp");
    std::cout << input << std::endl;
    double version = input["version"];                      // version = 1.1
    double pi = input["math"]["pi"];                        // pi = 3.14
    const std::vector<int>& fibo = input["math"]["series"]["fibo"]; // fibo = 2, 3, 5, 7, 11, 13, 17, 19
    for (const auto& f : fibo)
        std::cout << f << ", ";
    std::cout << std::endl;
    std::vector<int> binary = input["binary"];              // binary = 0, 1, 0, 1, 1, 1, 1, 0
    for (const auto& b : binary)
        std::cout << b << ", ";
    std::cout << std::endl;
    std::cout << version << ", " << pi << std::endl;
    std::vector<std::vector<int>> binary2D = input["binary"].matrix<int>();              // binary = 0, 1, 0, 1, 1, 1, 1, 0
    const std::vector<std::string> dir = input["dir"];
    std::cout << dir[0] << std::endl;
}