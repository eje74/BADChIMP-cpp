#include "LBsnippets.h"
#include <iomanip>
#include <sstream>
#include <iostream>
#ifdef __linux__ 
    #include <unistd.h>
#endif


void printAsciiToScreen(int nX, int nY, ScalarField& val, int** labels)
{
    unsigned int tmp;
    std::string gs = " .:-=+*#%@@@@@";

    for (int y = nY-1; y >= 0 ; y -= 2) {
        for (int x = 0; x < nX; ++x) {

            tmp = static_cast<unsigned int>(val(0, labels[y][x]) * 10);
            std::cout << std::setw(1) << gs[tmp];
        }
        // std::cout << std::endl << std::endl;
        std::cout << std::endl;
    }
}


void printAsciiToScreen(int nX, int nY, ScalarField& val, int** labels, double pauseInSeconds)
{
    //std::cout << "SCREEN" << std::endl << std::endl;
    unsigned int tmp;
    std::string gs = " .:-=+*#%@@@@@";
    unsigned int sleepMicroSec = static_cast<unsigned int>(1e6 * pauseInSeconds);
    for (int y = nY-1; y >= 0 ; y -= 2) {
        for (int x = 0; x < nX; ++x) {

            tmp = static_cast<unsigned int>(val(0, labels[y][x]) * 10);
            //std::cout << std::setw(5) << labels[y][x];
            std::cout << std::setw(1) << gs[tmp];
        }
    }
    std::cout << std::endl;
#ifdef __linux__ 
    usleep(sleepMicroSec);
#endif
    
}

