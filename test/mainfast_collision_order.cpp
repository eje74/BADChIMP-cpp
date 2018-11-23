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
#define N_ITERATIONS 10000
#define NX 250
#define NXQ 2250
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
#define DNY (NXQ + 18)  // pos = nx + ny*DNY
#define C2_INV 3.0
#define C2_INV_2 1.5
#define C4_INV 9.0
#define C4_INV_2 4.5

#define NEIG1 9
#define NEIG2 DNY
#define NEIG3 -9
#define NEIG4 -DNY

#define NEIG5 (9 + DNY)
#define NEIG6 (DNY - 9)
#define NEIG7 (-9 - DNY)
#define NEIG8 (9 - DNY)

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
    int pos_rho;
    double uu, cu, uF;

    int FIELD_SIZE = (NX + 2) * (NY + 2);

    double * data = new double [21 * FIELD_SIZE];

    double* F_EVEN = data;
    double* F_ODD = data + 9 * FIELD_SIZE;

    double* RHO = data + 18 * FIELD_SIZE;
    double* VX = data + 19 * FIELD_SIZE;
    double* VY = data + 20 * FIELD_SIZE;



    // INITIATE FIELDS
    pos = 9 + DNY;
    pos_rho = 3 + NX; 
    for (int j = 1; j <= NY; j++) {
        for (int i = 1; i <= NX; i++) {
            // Macroscopic
            RHO[pos_rho] = 1.0;
            VX[pos_rho] = 0.0;
            VY[pos_rho] = 0.0;
            // Microscopic
            F_EVEN[pos] = W0 * RHO[pos_rho];

            F_EVEN[pos + 1] = W1 * ( RHO[pos_rho] - 0.5 * C2_INV * CF1 );
            F_EVEN[pos + 2] = W1 * ( RHO[pos_rho] - 0.5 * C2_INV * CF2 );
            F_EVEN[pos + 3] = W1 * ( RHO[pos_rho] - 0.5 * C2_INV * CF3 );
            F_EVEN[pos + 4] = W1 * ( RHO[pos_rho] - 0.5 * C2_INV * CF4 );

            F_EVEN[pos + 5] = W5 * ( RHO[pos_rho] - 0.5 * C2_INV * CF5 );
            F_EVEN[pos + 6] = W5 * ( RHO[pos_rho] - 0.5 * C2_INV * CF6 );
            F_EVEN[pos + 7] = W5 * ( RHO[pos_rho] - 0.5 * C2_INV * CF7 );
            F_EVEN[pos + 8] = W5 * ( RHO[pos_rho] - 0.5 * C2_INV * CF8 );

            pos += 9;
	    pos_rho += 1;
        }
        pos += 18;
	pos_rho += 2;
    }

    // MAIN LOOP
    for (int n = 0; n < N_ITERATIONS; n += 2) {
        // EVEN FUNCTION
        pos = 9 + DNY;
	pos_rho = 3 + NX;
        for (int j = 1; j <= NY; j++) {
            for (int i = 1; i <= NX; i++) {
                // Rho + velocity
                RHO[pos_rho] = F_EVEN[pos] + F_EVEN[pos + 1] + F_EVEN[pos + 2] + F_EVEN[pos + 3] + F_EVEN[pos + 4] +
                                          F_EVEN[pos + 5] + F_EVEN[pos + 6] + F_EVEN[pos + 7] + F_EVEN[pos + 8];

                VX[pos_rho]  = (F_EVEN[pos + 1] - F_EVEN[pos + 3] +
                            F_EVEN[pos + 5] - F_EVEN[pos + 6] - F_EVEN[pos + 7] + F_EVEN[pos + 8] + 0.5 * FX) / RHO[pos_rho];

                VY[pos_rho]  = (F_EVEN[pos + 2] - F_EVEN[pos + 4] +
                            F_EVEN[pos + 5] + F_EVEN[pos + 6] - F_EVEN[pos + 7] - F_EVEN[pos + 8] + 0.5 * FY) / RHO[pos_rho];

                // Collision and propagation
                uu = VX[pos_rho] * VX[pos_rho] + VY[pos_rho] * VY[pos_rho];
                uF = VX[pos_rho] * FX + VY[pos_rho] * FY;

                F_ODD[pos] = (1.0 - OMEGA) * F_EVEN[pos]
                        + OMEGA * W0 * RHO[pos_rho] * (1.0 - C2_INV_2 * uu)
                        - (1.0 - 0.5*OMEGA) * W0 * C2_INV * uF;

                cu = VX[pos_rho];
                F_ODD[pos + 1 + NEIG1] = (1.0 - OMEGA) * F_EVEN[pos + 1]
                        + OMEGA * W1 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF1 + C4_INV * CF1 * cu - C2_INV * uF );

                cu = VY[pos_rho];
                F_ODD[pos + 2 + NEIG2] = (1.0 - OMEGA) * F_EVEN[pos + 2]
                        + OMEGA * W1 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF2 + C4_INV * CF2 * cu - C2_INV * uF );

                cu = -VX[pos_rho];
                F_ODD[pos + 3 + NEIG3] = (1.0 - OMEGA) * F_EVEN[pos + 3]
                        + OMEGA * W1 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF3 + C4_INV * CF3 * cu - C2_INV * uF );

                cu = -VY[pos_rho];
                F_ODD[pos + 4 + NEIG4] = (1.0 - OMEGA) * F_EVEN[pos + 4]
                        + OMEGA * W1 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF4 + C4_INV * CF4 * cu - C2_INV * uF );

                cu = VX[pos_rho] + VY[pos_rho];
                F_ODD[pos + 5 + NEIG5] = (1.0 - OMEGA) * F_EVEN[pos + 5]
                        + OMEGA * W5 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF5 + C4_INV * CF5 * cu - C2_INV * uF );

                cu = -VX[pos_rho] + VY[pos_rho];
                F_ODD[pos + 6 + NEIG6] = (1.0 - OMEGA) * F_EVEN[pos + 6]
                        + OMEGA * W5 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF6 + C4_INV * CF6 * cu - C2_INV * uF );

                cu = -VX[pos_rho] - VY[pos_rho];
                F_ODD[pos + 7 + NEIG7] = (1.0 - OMEGA) * F_EVEN[pos + 7]
                        + OMEGA * W5 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF7 + C4_INV * CF7 * cu - C2_INV * uF );

                cu =  VX[pos_rho] - VY[pos_rho];
                F_ODD[pos + 8 + NEIG8] = (1.0 - OMEGA) * F_EVEN[pos + 8]
                        + OMEGA * W5 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF8 + C4_INV * CF8 * cu - C2_INV * uF );

                pos += 9;
		pos_rho += 1;
	    } // END nx
            pos += 18;
	    pos_rho += 2;
        } // END ny

        // EVEN BOUNDARY CONDITIONS
        // -- Periodic left boundary
        pos = 9 + DNY;
	pos_rho = 3 + NX;
        F_ODD[pos + 1] = F_ODD[pos + 1 + NXQ];
        F_ODD[pos + 8] = F_ODD[pos + 8 + NXQ];

        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F_ODD[pos + 1] = F_ODD[pos + 1 + NXQ];
            F_ODD[pos + 5] = F_ODD[pos + 5 + NXQ];
            F_ODD[pos + 8] = F_ODD[pos + 8 + NXQ];
        }

        pos += DNY;
        F_ODD[pos + 1] = F_ODD[pos + 1 + NXQ];
        F_ODD[pos + 5] = F_ODD[pos + 5 + NXQ];

        // -- Periodic right boundary
        pos = NXQ + DNY;
        F_ODD[pos + 3] = F_ODD[pos + 3 - NXQ];
        F_ODD[pos + 7] = F_ODD[pos + 7 - NXQ];
        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F_ODD[pos + 3] = F_ODD[pos + 3 - NXQ];
            F_ODD[pos + 6] = F_ODD[pos + 6 - NXQ];
            F_ODD[pos + 7] = F_ODD[pos + 7 - NXQ];
        }
        pos += DNY;
        F_ODD[pos + 3] = F_ODD[pos + 3 - NXQ];
        F_ODD[pos + 6] = F_ODD[pos + 6 - NXQ];

        // -- Wall bottom
        pos = DNY;
        for (int i = 1; i <= NX; i++) {
            pos += 9;
            F_ODD[pos + 2] = F_ODD[pos + 4 + NEIG4];
            F_ODD[pos + 5] = F_ODD[pos + 7 + NEIG7];
            F_ODD[pos + 6] = F_ODD[pos + 8 + NEIG8];
        }
        // -- Wall top
        pos = DNY * NY;
        for (int i = 1; i <= NX; i++) {
            pos += 9;
            F_ODD[pos + 4] = F_ODD[pos + 2 + NEIG2];
            F_ODD[pos + 7] = F_ODD[pos + 5 + NEIG5];
            F_ODD[pos + 8] = F_ODD[pos + 6 + NEIG6];
        }

        // ODD FUNCTIONS
        pos = 9 + DNY;
	pos_rho = 3 + NX;
        for (int j = 1; j <= NY; j++) {
            for (int i = 1; i <= NX; i++) {
                // Rho + velocity
                RHO[pos_rho] = F_ODD[pos] + F_ODD[pos + 1] + F_ODD[pos + 2] + F_ODD[pos + 3] + F_ODD[pos + 4] +
                                          F_ODD[pos + 5] + F_ODD[pos + 6] + F_ODD[pos + 7] + F_ODD[pos + 8];

                VX[pos_rho]  = (F_ODD[pos + 1] - F_ODD[pos + 3] +
                            F_ODD[pos + 5] - F_ODD[pos + 6] - F_ODD[pos + 7] + F_ODD[pos + 8] + 0.5 * FX) / RHO[pos_rho];

                VY[pos_rho]  = (F_ODD[pos + 2] - F_ODD[pos + 4] +
                            F_ODD[pos + 5] + F_ODD[pos + 6] - F_ODD[pos +7] - F_ODD[pos + 8] + 0.5 * FY) / RHO[pos_rho];

                // Collision and propagation
                uu = VX[pos_rho] * VX[pos_rho] + VY[pos_rho] * VY[pos_rho];
                uF = VX[pos_rho] * FX + VY[pos_rho] * FY;

                F_EVEN[pos] = (1.0 - OMEGA) * F_ODD[pos]
                        + OMEGA * W0 * RHO[pos_rho] * (1.0 - C2_INV_2 * uu)
                        - (1.0 - 0.5*OMEGA) * W0 * C2_INV * uF;

                cu = VX[pos_rho];
                F_EVEN[pos + 1 + NEIG1] = (1.0 - OMEGA) * F_ODD[pos + 1]
                        + OMEGA * W1 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF1 + C4_INV * CF1 * cu - C2_INV * uF );

                cu = VY[pos_rho];
                F_EVEN[pos + 2 + NEIG2] = (1.0 - OMEGA) * F_ODD[pos + 2]
                        + OMEGA * W1 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF2 + C4_INV * CF2 * cu - C2_INV * uF );

                cu = -VX[pos_rho];
                F_EVEN[pos  + 3 + NEIG3] = (1.0 - OMEGA) * F_ODD[pos + 3]
                        + OMEGA * W1 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF3 + C4_INV * CF3 * cu - C2_INV * uF );

                cu = -VY[pos_rho];
                F_EVEN[pos + 4 + NEIG4] = (1.0 - OMEGA) * F_ODD[pos + 4]
                        + OMEGA * W1 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF4 + C4_INV * CF4 * cu - C2_INV * uF );

                cu = VX[pos_rho] + VY[pos_rho];
                F_EVEN[pos + 5 + NEIG5] = (1.0 - OMEGA) * F_ODD[pos + 5]
                        + OMEGA * W5 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF5 + C4_INV * CF5 * cu - C2_INV * uF );

                cu = -VX[pos_rho] + VY[pos_rho];
                F_EVEN[pos + 6 + NEIG6] = (1.0 - OMEGA) * F_ODD[pos + 6]
                        + OMEGA * W5 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF6 + C4_INV * CF6 * cu - C2_INV * uF );

                cu = -VX[pos_rho] - VY[pos_rho];
                F_EVEN[pos + 7 + NEIG7] = (1.0 - OMEGA) * F_ODD[pos + 7]
                        + OMEGA * W5 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF7 + C4_INV * CF7 * cu - C2_INV * uF );

                cu =  VX[pos_rho] - VY[pos_rho];
                F_EVEN[pos + 8 + NEIG8] = (1.0 - OMEGA) * F_ODD[pos + 8]
                        + OMEGA * W5 * RHO[pos_rho] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF8 + C4_INV * CF8 * cu - C2_INV * uF );

                pos += 9;
		pos_rho += 1;
            } // END nx
            pos += 18;
	    pos_rho += 2;
        } // END ny

        // EVEN BOUNDARY CONDITIONS
        // -- Periodic left boundary
        pos = 9 + DNY;
        F_EVEN[pos + 1] = F_EVEN[pos + 1 + NXQ];
        F_EVEN[pos + 8] = F_EVEN[pos + 8 + NXQ];

        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F_EVEN[pos + 1] = F_EVEN[pos + 1 + NXQ];
            F_EVEN[pos + 5] = F_EVEN[pos + 5 + NXQ];
            F_EVEN[pos + 8] = F_EVEN[pos + 8 + NXQ];
        }

        pos += DNY;
        F_EVEN[pos + 1] = F_EVEN[pos + 1 + NXQ];
        F_EVEN[pos + 5] = F_EVEN[pos + 5 + NXQ];

        // -- Periodic right boundary
        pos = NXQ + DNY;
        F_EVEN[pos + 3] = F_EVEN[pos + 3 - NXQ];
        F_EVEN[pos + 7] = F_EVEN[pos + 7 - NXQ];
        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F_EVEN[pos + 3] = F_EVEN[pos + 3 - NXQ];
            F_EVEN[pos + 6] = F_EVEN[pos + 6 - NXQ];
            F_EVEN[pos + 7] = F_EVEN[pos + 7 - NXQ];
        }
        pos += DNY;
        F_EVEN[pos + 3] = F_EVEN[pos + 3 - NXQ];
        F_EVEN[pos + 6] = F_EVEN[pos + 6 - NXQ];

        // -- Wall bottom
        pos = DNY;
        for (int i = 1; i <= NX; i++) {
            pos += 9;
            F_EVEN[pos + 2] = F_EVEN[pos + 4 + NEIG4];
            F_EVEN[pos + 5] = F_EVEN[pos + 7 + NEIG7];
            F_EVEN[pos + 6] = F_EVEN[pos + 8 + NEIG8];
        }
        // -- Wall top
        pos = DNY * NY;
        for (int i = 1; i <= NX; i++) {
            pos += 9;
            F_EVEN[pos + 4] = F_EVEN[pos + 2 + NEIG2];
            F_EVEN[pos + 7] = F_EVEN[pos + 5 + NEIG5];
            F_EVEN[pos + 8] = F_EVEN[pos + 6 + NEIG6];
        }


        // PRINT TO FILE

    } // END ITER

    // Print so that the simultion loop will not optimized away...
    double sum = 0;
    pos = 9 + DNY;
    for (int j = 1; j <= NY; j++) {
        for (int i = 1; i <= NX; i++) {
	  sum += F_EVEN[pos + 7]; 
	    pos += 9;
	  
	} // END nx
	pos += 18;
    } // END ny
    
    std::cout << sum << std::endl;

    return 0;
}
