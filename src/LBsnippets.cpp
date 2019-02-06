#include "LBsnippets.h"
#include <iomanip>
#include <sstream>
#include <iostream>
#include <unistd.h>

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
    unsigned int tmp;
    std::string gs = " .:-=+*#%@@@@@";
    unsigned int sleepMicroSec = static_cast<unsigned int>(1e6 * pauseInSeconds);

    for (int y = nY-1; y >= 0 ; y -= 2) {
        for (int x = 0; x < nX; ++x) {

            tmp = static_cast<unsigned int>(val(0, labels[y][x]) * 10);
            std::cout << std::setw(1) << gs[tmp];


        }
        // std::cout << std::endl << std::endl;
        std::cout << std::endl;
    }

    usleep(sleepMicroSec);
}
