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
    double tmp2=0.0, tmp3=0.0;
    std::string gs = " .:-=+*#%@@@@@";
    unsigned int sleepMicroSec = static_cast<unsigned int>(1e6 * pauseInSeconds);
    //for (int y = nY-1; y >= 0 ; y -= 2) {
    for (int y = nY-1; y >= 0 ; y --) {
        for (int x = 0; x < nX; ++x) {

	  //tmp = static_cast<unsigned int>(val(0, labels[y][x]) * 10);
	  //std::cout << std::setw(1) << gs[tmp];
	    tmp2+=val(0, labels[y][x]);
	    tmp3+=val(1, labels[y][x]);

        }
        // std::cout << std::endl << std::endl;
        //std::cout << std::endl;
    }
    std::cout << "rho0Tot = " << tmp2 << std::endl;
    std::cout << "rho1Tot = " << tmp3 << std::endl;
    std::cout << "rhoTot = " << tmp2 +tmp3 << std::endl;
    std::cout << std::endl;
    usleep(sleepMicroSec);
}

void printAsciiToScreen(int nX, int nY, int nZ, int itr, ScalarField& val, int*** labels, double pauseInSeconds)
{
    unsigned int tmp;
    double tmp2=0.0, tmp3=0.0;
    std::string gs = " .:-=+*#%@@@@@";
    unsigned int sleepMicroSec = static_cast<unsigned int>(1e6 * pauseInSeconds);
    std::cout << "Time step = "<< itr << std::endl;
    for (int z = 0; z < nZ ; z ++) {
        std::cout << "z=" << z << std::endl;
        for (int y = nY-1; y >= 0 ; y--) {
            for (int x = 0; x < nX; ++x) {

                //tmp = static_cast<unsigned int>(val(0, labels[z][y][x]) * 10);
                //std::cout << std::setw(1) << gs[tmp];
                // std::cout << std::setw(5) << labels[z][y][x];
                tmp2+=val(0, labels[z][y][x]);

            }
            // std::cout << std::endl << std::endl;
            // std::cout << std::endl;
        }
        // std::cout << std::endl;
    }

    std::cout << "rho0Tot = " << tmp2 << std::endl;
    std::cout << std::endl;
    usleep(sleepMicroSec);
}
