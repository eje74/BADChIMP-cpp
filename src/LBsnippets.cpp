#include "LBsnippets.h"
#include <iomanip>
#include <sstream>
#include <iostream>
#include <unistd.h>


void setConstPressure(const lbBase_t rho0, const lbBase_t rho1, const lbBase_t rhoConst)
{
    lbBase_t q0, q1;

    lbBase_t rhoTot = rho0 + rho1;
    lbBase_t qTmp = 2*(rhoConst - rhoTot);

    if (qTmp > 0) { // Only add oil
    } else { // Remove both oil and
    }



}


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
    usleep(sleepMicroSec);
}

/* void printAsciiToScreen(int nX, int nY, int nZ, int itr, ScalarField& val, int*** labels, double pauseInSeconds)
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
} */
