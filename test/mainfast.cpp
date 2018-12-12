#include <iostream>


/*
 *  Basis directory
 *
 *      6   2   5
 *        \ | /
 *      3 - 0 - 1
 *        / | \
 *      7   4   8
 *
 *
 *  Compile with :
 *                 g++ -std=c++11 -O3 mainfast.cpp
 *
 *  A Simple baseline test code:
 *  2D Poiseuille flow with gravity
 */


// CONSTANTS
#define N_ITERATIONS 10000
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

#define NEIG0 0
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

    //    double** F_EVEN = new double* [9];
    // double** F_ODD = new double* [9];

    //  for (int q = 0; q < 9; ++q) {
	// F_EVEN[q] = data + q * FIELD_SIZE;
    //  F_ODD[q] = data + (q + 9) * FIELD_SIZE;
    //	}

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

    double* F_EVEN = F0_EVEN;
    double* F_ODD = F0_ODD; 
    double cul[9];
    double cfl[9];

    double w[] = {W0, W1, W1, W1, W1, W5, W5, W5, W5};
    int neig[] = {0, NEIG1, NEIG2, NEIG3, NEIG4, NEIG5, NEIG6, NEIG7, NEIG8};


    // INITIATE FIELDS
    pos = 1 + DNY;
    for (int j = 1; j <= NY; j++) {
        for (int i = 1; i <= NX; i++) {
            // Macroscopic
            RHO[pos] = 1.0;
            VX[pos] = 0.0;
            VY[pos] = 0.0;
            // Microscopic

            F_EVEN[pos] = W0 * RHO[pos];

            F_EVEN[pos + FIELD_SIZE] = W1 * ( RHO[pos] - 0.5 * C2_INV * CF1 );
            F_EVEN[pos + 2 * FIELD_SIZE] = W1 * ( RHO[pos] - 0.5 * C2_INV * CF2 );
            F_EVEN[pos + 3 * FIELD_SIZE] = W1 * ( RHO[pos] - 0.5 * C2_INV * CF3 );
            F_EVEN[pos + 4 * FIELD_SIZE] = W1 * ( RHO[pos] - 0.5 * C2_INV * CF4 );

            F_EVEN[pos + 5 * FIELD_SIZE] = W5 * ( RHO[pos] - 0.5 * C2_INV * CF5 );
            F_EVEN[pos + 6 * FIELD_SIZE] = W5 * ( RHO[pos] - 0.5 * C2_INV * CF6 );
            F_EVEN[pos + 7 * FIELD_SIZE] = W5 * ( RHO[pos] - 0.5 * C2_INV * CF7 );
            F_EVEN[pos + 8 * FIELD_SIZE] = W5 * ( RHO[pos] - 0.5 * C2_INV * CF8 );

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
                RHO[pos] = F_EVEN[pos] + F_EVEN[pos + FIELD_SIZE] + F_EVEN[pos + 2 * FIELD_SIZE] + F_EVEN[pos + 3 * FIELD_SIZE] + F_EVEN[pos + 4 * FIELD_SIZE] +
                                          F_EVEN[pos + 5 * FIELD_SIZE] + F_EVEN[pos + 6 * FIELD_SIZE] + F_EVEN[pos + 7 * FIELD_SIZE] + F_EVEN[pos + 8 * FIELD_SIZE];

                VX[pos]  = (F_EVEN[pos + FIELD_SIZE] - F_EVEN[pos + 3 * FIELD_SIZE] +
                            F_EVEN[pos + 5 * FIELD_SIZE] - F_EVEN[pos + 6 * FIELD_SIZE] - F_EVEN[pos + 7 * FIELD_SIZE] + F_EVEN[pos + 8 * FIELD_SIZE] + 0.5 * FX) / RHO[pos];

                VY[pos]  = (F_EVEN[pos + 2 * FIELD_SIZE] - F_EVEN[pos + 4 * FIELD_SIZE] +
                            F_EVEN[pos + 5 * FIELD_SIZE] + F_EVEN[pos + 6 * FIELD_SIZE] - F_EVEN[pos + 7 * FIELD_SIZE] - F_EVEN[pos + 8 * FIELD_SIZE] + 0.5 * FY) / RHO[pos];


                // Collision and propagation
                uu = VX[pos] * VX[pos] + VY[pos] * VY[pos];
                uF = VX[pos] * FX + VY[pos] * FY;

                F_ODD[pos] = (1.0 - OMEGA) * F_EVEN[pos]
                        + OMEGA * W0 * RHO[pos] * (1.0 - C2_INV_2 * uu)
                        - (1.0 - 0.5*OMEGA) * W0 * C2_INV * uF;

                cu = VX[pos];
                F_ODD[pos + FIELD_SIZE + NEIG1] = (1.0 - OMEGA) * F_EVEN[pos + FIELD_SIZE]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF1 + C4_INV * CF1 * cu - C2_INV * uF );

                cu = VY[pos];
                F_ODD[pos + 2 * FIELD_SIZE + NEIG2] = (1.0 - OMEGA) * F_EVEN[pos + 2* FIELD_SIZE]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF2 + C4_INV * CF2 * cu - C2_INV * uF );

                cu = -VX[pos];
                F_ODD[pos + 3 * FIELD_SIZE + NEIG3] = (1.0 - OMEGA) * F_EVEN[pos + 3 * FIELD_SIZE]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF3 + C4_INV * CF3 * cu - C2_INV * uF );

                cu = -VY[pos];
                F_ODD[pos + 4 * FIELD_SIZE + NEIG4] = (1.0 - OMEGA) * F_EVEN[pos + 4 * FIELD_SIZE]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF4 + C4_INV * CF4 * cu - C2_INV * uF );

                cu = VX[pos] + VY[pos];
                F_ODD[pos + 5 * FIELD_SIZE + NEIG5] = (1.0 - OMEGA) * F_EVEN[pos + 5 * FIELD_SIZE]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF5 + C4_INV * CF5 * cu - C2_INV * uF );

                cu = -VX[pos] + VY[pos];
                F_ODD[pos + 6 * FIELD_SIZE + NEIG6] = (1.0 - OMEGA) * F_EVEN[pos + 6 * FIELD_SIZE]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF6 + C4_INV * CF6 * cu - C2_INV * uF );

                cu = -VX[pos] - VY[pos];
                F_ODD[pos + 7 * FIELD_SIZE + NEIG7] = (1.0 - OMEGA) * F_EVEN[pos + 7 * FIELD_SIZE]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF7 + C4_INV * CF7 * cu - C2_INV * uF );

                cu =  VX[pos] - VY[pos];
                F_ODD[pos + 8 * FIELD_SIZE + NEIG8] = (1.0 - OMEGA) * F_EVEN[pos + 8 * FIELD_SIZE]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF8 + C4_INV * CF8 * cu - C2_INV * uF );

		/*		cul[0] = 0.0;
		cul[1] = VX[pos];
		cul[2] = VY[pos];
		cul[3] = -VX[pos];
		cul[4] = -VY[pos];
		cul[5] = VX[pos] + VY[pos];
		cul[6] = -VX[pos] + VY[pos];
		cul[7] = -VX[pos] - VY[pos];
		cul[8] =  VX[pos] - VY[pos];
		
		cfl[0] = 0.0;
		cfl[1] =  FX;
		cfl[2] =  FY;
		cfl[3] = -FX;
		cfl[4] = -FY;

		cfl[5] = ( FX + FY);
		cfl[6] = ( FY - FX);
		cfl[7] = (-FX - FY);
		cfl[8] = ( FX - FY);
		
		
		F_ODD[pos + NEIG0 + 0*FIELD_SIZE] = 
		     (1.0 - OMEGA) * F_EVEN[pos + 0*FIELD_SIZE]
		      + OMEGA * w[0] * RHO[pos] * (1.0 + C2_INV * cul[0] + C4_INV_2 * cul[0] * cul[0] - C2_INV_2 * uu)
		      + (1.0 - 0.5*OMEGA) * w[0] * (C2_INV * cfl[0] + C4_INV * cfl[0] * cul[0] - C2_INV * uF );
		F_ODD[pos + NEIG1 + 1*FIELD_SIZE] = 
		     (1.0 - OMEGA) * F_EVEN[pos + 1*FIELD_SIZE]
		      + OMEGA * w[1] * RHO[pos] * (1.0 + C2_INV * cul[1] + C4_INV_2 * cul[1] * cul[1] - C2_INV_2 * uu)
		      + (1.0 - 0.5*OMEGA) * w[1] * (C2_INV * cfl[1] + C4_INV * cfl[1] * cul[1] - C2_INV * uF );
		F_ODD[pos + NEIG2 + 2*FIELD_SIZE] = 
		     (1.0 - OMEGA) * F_EVEN[pos + 2*FIELD_SIZE]
		      + OMEGA * w[2] * RHO[pos] * (1.0 + C2_INV * cul[2] + C4_INV_2 * cul[2] * cul[2] - C2_INV_2 * uu)
		      + (1.0 - 0.5*OMEGA) * w[2] * (C2_INV * cfl[2] + C4_INV * cfl[2] * cul[2] - C2_INV * uF );
		F_ODD[pos + NEIG3 + 3*FIELD_SIZE] = 
		     (1.0 - OMEGA) * F_EVEN[pos + 3*FIELD_SIZE]
		      + OMEGA * w[3] * RHO[pos] * (1.0 + C2_INV * cul[3] + C4_INV_2 * cul[3] * cul[3] - C2_INV_2 * uu)
		      + (1.0 - 0.5*OMEGA) * w[3] * (C2_INV * cfl[3] + C4_INV * cfl[3] * cul[3] - C2_INV * uF );
		F_ODD[pos + NEIG4 + 4*FIELD_SIZE] = 
		     (1.0 - OMEGA) * F_EVEN[pos + 4*FIELD_SIZE]
		      + OMEGA * w[4] * RHO[pos] * (1.0 + C2_INV * cul[4] + C4_INV_2 * cul[4] * cul[4] - C2_INV_2 * uu)
		      + (1.0 - 0.5*OMEGA) * w[4] * (C2_INV * cfl[4] + C4_INV * cfl[4] * cul[4] - C2_INV * uF );
		F_ODD[pos + NEIG5 + 5*FIELD_SIZE] = 
		     (1.0 - OMEGA) * F_EVEN[pos + 5*FIELD_SIZE]
		      + OMEGA * w[5] * RHO[pos] * (1.0 + C2_INV * cul[5] + C4_INV_2 * cul[5] * cul[5] - C2_INV_2 * uu)
		      + (1.0 - 0.5*OMEGA) * w[5] * (C2_INV * cfl[5] + C4_INV * cfl[5] * cul[5] - C2_INV * uF );
		F_ODD[pos + NEIG6 + 6*FIELD_SIZE] = 
		     (1.0 - OMEGA) * F_EVEN[pos + 6*FIELD_SIZE]
		      + OMEGA * w[6] * RHO[pos] * (1.0 + C2_INV * cul[6] + C4_INV_2 * cul[6] * cul[6] - C2_INV_2 * uu)
		      + (1.0 - 0.5*OMEGA) * w[6] * (C2_INV * cfl[6] + C4_INV * cfl[6] * cul[6] - C2_INV * uF );
		F_ODD[pos + NEIG7 + 7*FIELD_SIZE] = 
		     (1.0 - OMEGA) * F_EVEN[pos + 7*FIELD_SIZE]
		      + OMEGA * w[7] * RHO[pos] * (1.0 + C2_INV * cul[7] + C4_INV_2 * cul[7] * cul[7] - C2_INV_2 * uu)
		      + (1.0 - 0.5*OMEGA) * w[7] * (C2_INV * cfl[7] + C4_INV * cfl[7] * cul[7] - C2_INV * uF );
		F_ODD[pos + NEIG8 + 8*FIELD_SIZE] = 
		     (1.0 - OMEGA) * F_EVEN[pos + 8*FIELD_SIZE]
		      + OMEGA * w[8] * RHO[pos] * (1.0 + C2_INV * cul[8] + C4_INV_2 * cul[8] * cul[8] - C2_INV_2 * uu)
		      + (1.0 - 0.5*OMEGA) * w[8] * (C2_INV * cfl[8] + C4_INV * cfl[8] * cul[8] - C2_INV * uF ); */
		
                pos += 1;
            } // END nx
            pos += 2;
        } // END ny

        // EVEN BOUNDARY CONDITIONS
        // -- Periodic left boundary
        pos = 1 + DNY;
        F_ODD[pos + FIELD_SIZE] = F_ODD[pos + FIELD_SIZE + NX];
        F_ODD[pos + 8 * FIELD_SIZE] = F_ODD[pos + 8 * FIELD_SIZE + NX];

        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F_ODD[pos + FIELD_SIZE] = F_ODD[pos + FIELD_SIZE + NX];
            F_ODD[pos + 5 * FIELD_SIZE] = F_ODD[pos + 5 * FIELD_SIZE + NX];
            F_ODD[pos + 8 * FIELD_SIZE] = F_ODD[pos + 8 * FIELD_SIZE + NX];
        }

        pos += DNY;
        F_ODD[pos + FIELD_SIZE] = F_ODD[pos + FIELD_SIZE + NX];
        F_ODD[pos + 5 * FIELD_SIZE] = F_ODD[pos + 5 * FIELD_SIZE + NX];

        // -- Periodic right boundary
        pos = NX + DNY;
        F_ODD[pos + 3 * FIELD_SIZE] = F_ODD[pos + 3 * FIELD_SIZE - NX];
        F_ODD[pos + 7 * FIELD_SIZE] = F_ODD[pos + 7 * FIELD_SIZE - NX];
        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F_ODD[pos + 3 * FIELD_SIZE] = F_ODD[pos + 3 * FIELD_SIZE - NX];
            F_ODD[pos + 6 * FIELD_SIZE] = F_ODD[pos + 6 * FIELD_SIZE - NX];
            F_ODD[pos + 7 * FIELD_SIZE] = F_ODD[pos + 7 * FIELD_SIZE - NX];
        }
        pos += DNY;
        F_ODD[pos + 3 * FIELD_SIZE] = F_ODD[pos + 3 * FIELD_SIZE - NX];
        F_ODD[pos + 6 * FIELD_SIZE] = F_ODD[pos + 6 * FIELD_SIZE - NX];

        // -- Wall bottom
        pos = DNY;
        for (int i = 1; i <= NX; i++) {
            pos += 1;
            F_ODD[pos + 2* FIELD_SIZE] = F_ODD[pos + 4 * FIELD_SIZE + NEIG4];
            F_ODD[pos + 5 * FIELD_SIZE] = F_ODD[pos + 7 * FIELD_SIZE + NEIG7];
            F_ODD[pos + 6 * FIELD_SIZE] = F_ODD[pos + 8 * FIELD_SIZE + NEIG8];
        }
        // -- Wall top
        pos = DNY * NY;
        for (int i = 1; i <= NX; i++) {
            pos += 1;
            F_ODD[pos + 4 * FIELD_SIZE] = F_ODD[pos + 2* FIELD_SIZE + NEIG2];
            F_ODD[pos + 7 * FIELD_SIZE] = F_ODD[pos + 5 * FIELD_SIZE + NEIG5];
            F_ODD[pos + 8 * FIELD_SIZE] = F_ODD[pos + 6 * FIELD_SIZE + NEIG6];
        }

        // ODD FUNCTIONS
        pos = 1 + DNY;
        for (int j = 1; j <= NY; j++) {
            for (int i = 1; i <= NX; i++) {
                // Rho + velocity
                RHO[pos] = F_ODD[pos] + F_ODD[pos + FIELD_SIZE] + F_ODD[pos + 2* FIELD_SIZE] + F_ODD[pos + 3 * FIELD_SIZE] + F_ODD[pos + 4 * FIELD_SIZE] +
                                          F_ODD[pos + 5 * FIELD_SIZE] + F_ODD[pos + 6 * FIELD_SIZE] + F_ODD[pos + 7 * FIELD_SIZE] + F_ODD[pos + 8 * FIELD_SIZE];

                VX[pos]  = (F_ODD[pos + FIELD_SIZE] - F_ODD[pos + 3 * FIELD_SIZE] +
                            F_ODD[pos + 5 * FIELD_SIZE] - F_ODD[pos + 6 * FIELD_SIZE] - F_ODD[pos + 7 * FIELD_SIZE] + F_ODD[pos + 8 * FIELD_SIZE] + 0.5 * FX) / RHO[pos];

                VY[pos]  = (F_ODD[pos + 2* FIELD_SIZE] - F_ODD[pos + 4 * FIELD_SIZE] +
                            F_ODD[pos + 5 * FIELD_SIZE] + F_ODD[pos + 6 * FIELD_SIZE] - F_ODD[pos + 7 * FIELD_SIZE] - F_ODD[pos + 8 * FIELD_SIZE] + 0.5 * FY) / RHO[pos];


                // Collision and propagation
                uu = VX[pos] * VX[pos] + VY[pos] * VY[pos];
                uF = VX[pos] * FX + VY[pos] * FY;

                F_EVEN[pos] = (1.0 - OMEGA) * F_ODD[pos]
                        + OMEGA * W0 * RHO[pos] * (1.0 - C2_INV_2 * uu)
                        - (1.0 - 0.5*OMEGA) * W0 * C2_INV * uF;

                cu = VX[pos];
                F_EVEN[pos + FIELD_SIZE + NEIG1] = (1.0 - OMEGA) * F_ODD[pos + FIELD_SIZE]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF1 + C4_INV * CF1 * cu - C2_INV * uF );

                cu = VY[pos];
                F_EVEN[pos + 2 * FIELD_SIZE + NEIG2] = (1.0 - OMEGA) * F_ODD[pos + 2* FIELD_SIZE]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF2 + C4_INV * CF2 * cu - C2_INV * uF );

                cu = -VX[pos];
                F_EVEN[pos + 3 * FIELD_SIZE + NEIG3] = (1.0 - OMEGA) * F_ODD[pos + 3 * FIELD_SIZE]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF3 + C4_INV * CF3 * cu - C2_INV * uF );

                cu = -VY[pos];
                F_EVEN[pos + 4 * FIELD_SIZE + NEIG4] = (1.0 - OMEGA) * F_ODD[pos + 4 * FIELD_SIZE]
                        + OMEGA * W1 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W1 * (C2_INV * CF4 + C4_INV * CF4 * cu - C2_INV * uF );

                cu = VX[pos] + VY[pos];
                F_EVEN[pos + 5 * FIELD_SIZE + NEIG5] = (1.0 - OMEGA) * F_ODD[pos + 5 * FIELD_SIZE]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF5 + C4_INV * CF5 * cu - C2_INV * uF );

                cu = -VX[pos] + VY[pos];
                F_EVEN[pos + 6 * FIELD_SIZE + NEIG6] = (1.0 - OMEGA) * F_ODD[pos + 6 * FIELD_SIZE]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF6 + C4_INV * CF6 * cu - C2_INV * uF );

                cu = -VX[pos] - VY[pos];
                F_EVEN[pos + 7 * FIELD_SIZE + NEIG7] = (1.0 - OMEGA) * F_ODD[pos + 7 * FIELD_SIZE]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF7 + C4_INV * CF7 * cu - C2_INV * uF );

                cu =  VX[pos] - VY[pos];
                F_EVEN[pos + 8 * FIELD_SIZE + NEIG8] = (1.0 - OMEGA) * F_ODD[pos + 8 * FIELD_SIZE]
                        + OMEGA * W5 * RHO[pos] * (1.0 + C2_INV * cu + C4_INV_2 * cu * cu - C2_INV_2 * uu)
                        + (1.0 - 0.5*OMEGA) * W5 * (C2_INV * CF8 + C4_INV * CF8 * cu - C2_INV * uF );

                pos += 1;
            } // END nx
            pos += 2;
        } // END ny

        // EVEN BOUNDARY CONDITIONS
        // -- Periodic left boundary
        pos = 1 + DNY;
        F_EVEN[pos + FIELD_SIZE] = F_EVEN[pos + FIELD_SIZE + NX];
        F_EVEN[pos + 8 * FIELD_SIZE] = F_EVEN[pos + 8 * FIELD_SIZE + NX];

        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F_EVEN[pos + FIELD_SIZE] = F_EVEN[pos + FIELD_SIZE + NX];
            F_EVEN[pos + 5 * FIELD_SIZE] = F_EVEN[pos + 5 * FIELD_SIZE + NX];
            F_EVEN[pos + 8 * FIELD_SIZE] = F_EVEN[pos + 8 * FIELD_SIZE + NX];
        }

        pos += DNY;
        F_EVEN[pos + FIELD_SIZE] = F_EVEN[pos + FIELD_SIZE + NX];
        F_EVEN[pos + 5 * FIELD_SIZE] = F_EVEN[pos + 5 * FIELD_SIZE + NX];

        // -- Periodic right boundary
        pos = NX + DNY;
        F_EVEN[pos + 3 * FIELD_SIZE] = F_EVEN[pos + 3 * FIELD_SIZE - NX];
        F_EVEN[pos + 7 * FIELD_SIZE] = F_EVEN[pos + 7 * FIELD_SIZE - NX];
        for (int j = 2; j < NY; j++) {
            pos += DNY;
            F_EVEN[pos + 3 * FIELD_SIZE] = F_EVEN[pos + 3 * FIELD_SIZE - NX];
            F_EVEN[pos + 6 * FIELD_SIZE] = F_EVEN[pos + 6 * FIELD_SIZE - NX];
            F_EVEN[pos + 7 * FIELD_SIZE] = F_EVEN[pos + 7 * FIELD_SIZE - NX];
        }
        pos += DNY;
        F_EVEN[pos + 3 * FIELD_SIZE] = F_EVEN[pos + 3 * FIELD_SIZE - NX];
        F_EVEN[pos + 6 * FIELD_SIZE] = F_EVEN[pos + 6 * FIELD_SIZE - NX];

        // -- Wall bottom
        pos = DNY;
        for (int i = 1; i <= NX; i++) {
            pos += 1;
            F_EVEN[pos + 2 * FIELD_SIZE] = F_EVEN[pos + 4 * FIELD_SIZE + NEIG4];
            F_EVEN[pos + 5 * FIELD_SIZE] = F_EVEN[pos + 7 * FIELD_SIZE + NEIG7];
            F_EVEN[pos + 6 * FIELD_SIZE] = F_EVEN[pos + 8 * FIELD_SIZE + NEIG8];
        }
        // -- Wall top
        pos = DNY * NY;
        for (int i = 1; i <= NX; i++) {
            pos += 1;
            F_EVEN[pos + 4 * FIELD_SIZE] = F_EVEN[pos + 2 * FIELD_SIZE + NEIG2];
            F_EVEN[pos + 7 * FIELD_SIZE] = F_EVEN[pos + 5 * FIELD_SIZE + NEIG5];
            F_EVEN[pos + 8 * FIELD_SIZE] = F_EVEN[pos + 6 * FIELD_SIZE + NEIG6];
        }


        // PRINT TO FILE

    } // END ITER

    // Print so that the simultion loop will not optimized away...
    pos = 1 + DNY;
    std::cout << VX[pos] << std::endl;

    return 0;
}
