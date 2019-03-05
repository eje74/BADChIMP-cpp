#include <iostream>


//
// Basis directory
//
//      6   2   5
//        \ | /
//      3 - 0 - 1
//        / | \
//      7   4   8
//
//
//  Compile with :
//                 g++ -std=c++11 -O3 mainfast.cpp
//
//  A Simple baseline test code:
//  2D Poiseuille flow with gravity
//


// CONSTANTS
#define N_ITERATIONS 1000
#define NX 250
#define NY 100
#define OMEGA 1.3
#define FX 1.0e-8
#define FY 0.0

#define NU ((1.0/OMEGA - 0.5)/3.0)

// LATTICE
// -- Weights
#define W0 0.444444444444  // Rest particle
#define W1 0.111111111111  // Nearest neighbor
#define W5 0.027777777778  // Next nearest neighbor
#define DNY (NX + 2)  // pos = nx + ny*DNY
#define C2_INV 3.0
#define C2_INV_2 1.5
#define C4_INV 9.0
#define C4_INV_2 4.5

#define NEIG1 1
#define NEIG2 DNY
#define NEIG3 -1
#define NEIG4 -DNY

#define NEIG5 (1 + DNY)
#define NEIG6 (DNY - 1)
#define NEIG7 (-1 - DNY)
#define NEIG8 (1 - DNY)

#define CF1  FX
#define CF2  FY
#define CF3 -FX
#define CF4 -FY

#define CF5 ( FX + FY)
#define CF6 ( FY - FX)
#define CF7 (-FX - FY)
#define CF8 ( FX - FY)

// FIELDS
//#define FIELD_SIZE ( (NX + 2) * (NY + 2) )


int main()
{
    int pos;
    double uu, cu, uF;

    int FIELD_SIZE = (NX + 2) * (NY + 2);

    double * data = new double [21 * FIELD_SIZE];

    double* F0_EVEN = data;
    double* F1_EVEN = data + 1 * FIELD_SIZE;
    double* F2_EVEN = data + 2 * FIELD_SIZE;
    double* F3_EVEN = data + 3 * FIELD_SIZE;
    double* F4_EVEN = data + 4 * FIELD_SIZE;
    double* F5_EVEN = data + 5 * FIELD_SIZE;
    double* F6_EVEN = data + 6 * FIELD_SIZE;
    double* F7_EVEN = data + 7 * FIELD_SIZE;
    double* F8_EVEN = data + 8 * FIELD_SIZE;
    double* F0_ODD = data + 9 * FIELD_SIZE;
    double* F1_ODD = data + 10 * FIELD_SIZE;
    double* F2_ODD = data + 11 * FIELD_SIZE;
    double* F3_ODD = data + 12 * FIELD_SIZE;
    double* F4_ODD = data + 13 * FIELD_SIZE;
    double* F5_ODD = data + 14 * FIELD_SIZE;
    double* F6_ODD = data + 15 * FIELD_SIZE;
    double* F7_ODD = data + 16 * FIELD_SIZE;
    double* F8_ODD = data + 17 * FIELD_SIZE;

    double* RHO = data + 18 * FIELD_SIZE;
    double* VX = data + 19 * FIELD_SIZE;
    double* VY = data + 20 * FIELD_SIZE;



    // INITIATE FIELDS
    pos = 1 + DNY;
    for (int j = 1; j <= NY; j++) {
        for (int i = 1; i <= NX; i++) {
            // Macroscopic
            RHO[pos] = 1.0;
            VX[pos] = 0.0;
            VY[pos] = 0.0;
            // Microscopic
            F0_EVEN[pos] = W0 * RHO[pos];

            F1_EVEN[pos] = W1 * ( RHO[pos] - 0.5 * C2_INV * CF1 );
            F2_EVEN[pos] = W1 * ( RHO[pos] - 0.5 * C2_INV * CF2 );
            F3_EVEN[pos] = W1 * ( RHO[pos] - 0.5 * C2_INV * CF3 );
            F4_EVEN[pos] = W1 * ( RHO[pos] - 0.5 * C2_INV * CF4 );

            F5_EVEN[pos] = W5 * ( RHO[pos] - 0.5 * C2_INV * CF5 );
            F6_EVEN[pos] = W5 * ( RHO[pos] - 0.5 * C2_INV * CF6 );
            F7_EVEN[pos] = W5 * ( RHO[pos] - 0.5 * C2_INV * CF7 );
            F8_EVEN[pos] = W5 * ( RHO[pos] - 0.5 * C2_INV * CF8 );

            pos += 1;

        }
        pos += 2;
    }

    // MAIN LOOP
    for (int n = 0; n < N_ITERATIONS; n += 2) {
        // EVEN FUNCTION
        pos = 1 + DNY;
        for (int j = 1; j <= NY; j++) {
            for (int i = 1; i <= NX; i++) {
                // Rho + velocity
                RHO[pos] = F0_EVEN[pos] + F1_EVEN[pos] + F2_EVEN[pos] + F3_EVEN[pos] + F4_EVEN[pos] +
                                          F5_EVEN[pos] + F6_EVEN[pos] + F7_EVEN[pos] + F8_EVEN[pos];

                VX[pos]  = (F1_EVEN[pos] - F3_EVEN[pos] +
                            F5_EVEN[pos] - F6_EVEN[pos] - F7_EVEN[pos] + F8_EVEN[pos] + 0.5 * FX) / RHO[pos];

                VY[pos]  = (F2_EVEN[pos] - F4_EVEN[pos] +
                            F5_EVEN[pos] + F6_EVEN[pos] - F7_EVEN[pos] - F8_EVEN[pos] + 0.5 * FY) / RHO[pos];

                // Collision and propagation
                uu = VX[pos] * VX[pos] + VY[pos] * VY[pos];
                uF = VX[pos] * FX + VY[pos] * FY;

                F0_ODD[pos] = (1.0 - OMEGA) * F0_EVEN[pos]
                        + OMEGA * W0 * RHO[pos] * (1.0 - C2_INV_2 * uu)
                        - (1.0 - 0.5*OMEGA) * W0 * C2_INV * uF;

                cu = VX[pos];
                F1_ODD[pos + NEIG1] = (1.0 - OMEGA) * F1_EVEN[pos]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF1 + C4_INV * CF1 * cu - C2_INV * uF );

                cu = VY[pos];
                F2_ODD[pos + NEIG2] = (1.0 - OMEGA) * F2_EVEN[pos]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF2 + C4_INV * CF2 * cu - C2_INV * uF );

                cu = -VX[pos];
                F3_ODD[pos + NEIG3] = (1.0 - OMEGA) * F3_EVEN[pos]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF3 + C4_INV * CF3 * cu - C2_INV * uF );

                cu = -VY[pos];
                F4_ODD[pos + NEIG4] = (1.0 - OMEGA) * F4_EVEN[pos]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF4 + C4_INV * CF4 * cu - C2_INV * uF );

                cu = VX[pos] + VY[pos];
                F5_ODD[pos + NEIG5] = (1.0 - OMEGA) * F5_EVEN[pos]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF5 + C4_INV * CF5 * cu - C2_INV * uF );

                cu = -VX[pos] + VY[pos];
                F6_ODD[pos + NEIG6] = (1.0 - OMEGA) * F6_EVEN[pos]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF6 + C4_INV * CF6 * cu - C2_INV * uF );

                cu = -VX[pos] - VY[pos];
                F7_ODD[pos + NEIG7] = (1.0 - OMEGA) * F7_EVEN[pos]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF7 + C4_INV * CF7 * cu - C2_INV * uF );

                cu =  VX[pos] - VY[pos];
                F8_ODD[pos + NEIG8] = (1.0 - OMEGA) * F8_EVEN[pos]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF8 + C4_INV * CF8 * cu - C2_INV * uF );

                pos += 1;
            } // END nx
            pos += 2;
        } // END ny

        // EVEN BOUNDARY CONDITIONS
        // -- Periodic left boundary
        pos = 1 + DNY;
        F1_ODD[pos] = F1_ODD[pos + NX];
        F8_ODD[pos] = F8_ODD[pos + NX];

        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F1_ODD[pos] = F1_ODD[pos + NX];
            F5_ODD[pos] = F5_ODD[pos + NX];
            F8_ODD[pos] = F8_ODD[pos + NX];
        }

        pos += DNY;
        F1_ODD[pos] = F1_ODD[pos + NX];
        F5_ODD[pos] = F5_ODD[pos + NX];

        // -- Periodic right boundary
        pos = NX + DNY;
        F3_ODD[pos] = F3_ODD[pos - NX];
        F7_ODD[pos] = F7_ODD[pos - NX];
        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F3_ODD[pos] = F3_ODD[pos - NX];
            F6_ODD[pos] = F6_ODD[pos - NX];
            F7_ODD[pos] = F7_ODD[pos - NX];
        }
        pos += DNY;
        F3_ODD[pos] = F3_ODD[pos - NX];
        F6_ODD[pos] = F6_ODD[pos - NX];

        // -- Wall bottom
        pos = DNY;
        for (int i = 1; i <= NX; i++) {
            pos += 1;
            F2_ODD[pos] = F4_ODD[pos + NEIG4];
            F5_ODD[pos] = F7_ODD[pos + NEIG7];
            F6_ODD[pos] = F8_ODD[pos + NEIG8];
        }
        // -- Wall top
        pos = DNY * NY;
        for (int i = 1; i <= NX; i++) {
            pos += 1;
            F4_ODD[pos] = F2_ODD[pos + NEIG2];
            F7_ODD[pos] = F5_ODD[pos + NEIG5];
            F8_ODD[pos] = F6_ODD[pos + NEIG6];
        }

        // ODD FUNCTIONS
        pos = 1 + DNY;
        for (int j = 1; j <= NY; j++) {
            for (int i = 1; i <= NX; i++) {
                // Rho + velocity
                RHO[pos] = F0_ODD[pos] + F1_ODD[pos] + F2_ODD[pos] + F3_ODD[pos] + F4_ODD[pos] +
                                          F5_ODD[pos] + F6_ODD[pos] + F7_ODD[pos] + F8_ODD[pos];

                VX[pos]  = (F1_ODD[pos] - F3_ODD[pos] +
                            F5_ODD[pos] - F6_ODD[pos] - F7_ODD[pos] + F8_ODD[pos] + 0.5 * FX) / RHO[pos];

                VY[pos]  = (F2_ODD[pos] - F4_ODD[pos] +
                            F5_ODD[pos] + F6_ODD[pos] - F7_ODD[pos] - F8_ODD[pos] + 0.5 * FY) / RHO[pos];

                // Collision and propagation
                uu = VX[pos] * VX[pos] + VY[pos] * VY[pos];
                uF = VX[pos] * FX + VY[pos] * FY;

                F0_EVEN[pos] = (1.0 - OMEGA) * F0_ODD[pos]
                        + OMEGA * W0 * RHO[pos] * (1.0 - C2_INV_2 * uu)
                        - (1.0 - 0.5*OMEGA) * W0 * C2_INV * uF;

                cu = VX[pos];
                F1_EVEN[pos + NEIG1] = (1.0 - OMEGA) * F1_ODD[pos]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF1 + C4_INV * CF1 * cu - C2_INV * uF );

                cu = VY[pos];
                F2_EVEN[pos + NEIG2] = (1.0 - OMEGA) * F2_ODD[pos]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF2 + C4_INV * CF2 * cu - C2_INV * uF );

                cu = -VX[pos];
                F3_EVEN[pos + NEIG3] = (1.0 - OMEGA) * F3_ODD[pos]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF3 + C4_INV * CF3 * cu - C2_INV * uF );

                cu = -VY[pos];
                F4_EVEN[pos + NEIG4] = (1.0 - OMEGA) * F4_ODD[pos]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF4 + C4_INV * CF4 * cu - C2_INV * uF );

                cu = VX[pos] + VY[pos];
                F5_EVEN[pos + NEIG5] = (1.0 - OMEGA) * F5_ODD[pos]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF5 + C4_INV * CF5 * cu - C2_INV * uF );

                cu = -VX[pos] + VY[pos];
                F6_EVEN[pos + NEIG6] = (1.0 - OMEGA) * F6_ODD[pos]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF6 + C4_INV * CF6 * cu - C2_INV * uF );

                cu = -VX[pos] - VY[pos];
                F7_EVEN[pos + NEIG7] = (1.0 - OMEGA) * F7_ODD[pos]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF7 + C4_INV * CF7 * cu - C2_INV * uF );

                cu =  VX[pos] - VY[pos];
                F8_EVEN[pos + NEIG8] = (1.0 - OMEGA) * F8_ODD[pos]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF8 + C4_INV * CF8 * cu - C2_INV * uF );

                pos += 1;
            } // END nx
            pos += 2;
        } // END ny


        // EVEN BOUNDARY CONDITIONS
        // -- Periodic left boundary
        pos = 1 + DNY;
        F1_EVEN[pos] = F1_EVEN[pos + NX];
        F8_EVEN[pos] = F8_EVEN[pos + NX];

        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F1_EVEN[pos] = F1_EVEN[pos + NX];
            F5_EVEN[pos] = F5_EVEN[pos + NX];
            F8_EVEN[pos] = F8_EVEN[pos + NX];
        }

        pos += DNY;
        F1_EVEN[pos] = F1_EVEN[pos + NX];
        F5_EVEN[pos] = F5_EVEN[pos + NX];

        // -- Periodic right boundary
        pos = NX + DNY;
        F3_EVEN[pos] = F3_EVEN[pos - NX];
        F7_EVEN[pos] = F7_EVEN[pos - NX];
        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F3_EVEN[pos] = F3_EVEN[pos - NX];
            F6_EVEN[pos] = F6_EVEN[pos - NX];
            F7_EVEN[pos] = F7_EVEN[pos - NX];
        }
        pos += DNY;
        F3_EVEN[pos] = F3_EVEN[pos - NX];
        F6_EVEN[pos] = F6_EVEN[pos - NX];

        // -- Wall bottom
        pos = DNY;
        for (int i = 1; i <= NX; i++) {
            pos += 1;
            F2_EVEN[pos] = F4_EVEN[pos + NEIG4];
            F5_EVEN[pos] = F7_EVEN[pos + NEIG7];
            F6_EVEN[pos] = F8_EVEN[pos + NEIG8];
        }
        // -- Wall top
        pos = DNY * NY;
        for (int i = 1; i <= NX; i++) {
            pos += 1;
            F4_EVEN[pos] = F2_EVEN[pos + NEIG2];
            F7_EVEN[pos] = F5_EVEN[pos + NEIG5];
            F8_EVEN[pos] = F6_EVEN[pos + NEIG6];
        }


        // PRINT TO FILE

    } // END ITER

    // Print so that the simultion loop will not optimized away...
    pos = 1 + DNY;
    std::cout << VX[pos] << std::endl;

    return 0;
}
