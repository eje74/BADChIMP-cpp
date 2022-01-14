#include "src/io/Input.h"

int main() 
{
    Input input("input/input.dat");
    int max = input["iterations"]["max"];
    //std::cout << input["iterations"]["max"] << std::endl;
    std::cout << max << std::endl;
}