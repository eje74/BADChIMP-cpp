// //////////////////////////////////////////////
//
// BADChIMP std_case
//
// For documentation see:
//    doc/documentation.pdf
// 
// //////////////////////////////////////////////

#include "../LBSOLVER.h"
#include "../IO.h"
#include "./LBco2help.h"

// SET THE LATTICE TYPE
#define LT D2Q9

int main()
{
    // *********
    // SETUP MPI
    // *********
    MPI_Init(NULL, NULL);
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // ********************************
    // SETUP THE INPUT AND OUTPUT PATHS
    // ********************************
    std::string chimpDir = "./";
    // std::string chimpDir = "/home/AD.NORCERESEARCH.NO/esje/Programs/GitHub/BADCHiMP/";
    std::string mpiDir = chimpDir + "input/mpi/";
    std::string inputDir = chimpDir + "input/";
    std::string outputDir = chimpDir + "output/";

    // ***********************
    // SETUP GRID AND GEOMETRY
    // ***********************
    Input input(inputDir + "input.dat");
    LBvtk<LT> vtklb(mpiDir + "tmp" + std::to_string(myRank) + ".vtklb");
    Grid<LT> grid(vtklb);
    Nodes<LT> nodes(vtklb, grid);
    BndMpi<LT> mpiBoundary(vtklb, nodes, grid);
    // Set bulk nodes
    std::vector<int> bulkNodes = findBulkNodes(nodes);
    // Set solid boundary 
    std::vector<int> solidBoundaryNodes = findSolidBndNodes(nodes);

    // ***************
    // READ FROM INPUT
    // ***************
    // Number of iterations
    // int nIterations = static_cast<int>( input["iterations"]["max"]);
    int nIterations = input["iterations"]["max"];
    // Write interval
    // int nItrWrite = static_cast<int>( input["iterations"]["write"]);
    int nItrWrite = input["iterations"]["write"];

    // Fluid Flow
    // **********
    // Number of fluid fields
    const int nFluidFields = input["fluid"]["numfluids"];
    std::cout << "nFluidFields = " << nFluidFields << std::endl;
    // Relaxation time
    // lbBase_t tau = 1;
    // Kinematic viscosities
    std::valarray<lbBase_t> kin_visc = input["fluid"]["viscosity"];
    std::valarray<lbBase_t> kin_visc_inv(nFluidFields);
    for (int n = 0; n < nFluidFields; ++n){
      kin_visc_inv[n] = 1./kin_visc[n];
    }
    // Body force
    //VectorField<LT> bodyForce(1, 1);
    //bodyForce.set(0, 0) = input["fluid"]["bodyforce"].valarray();
    VectorField<LT> bodyForce(1, 1, input["fluid"]["bodyforce"]);
    std::cout << "bodyForce\n";
    for (const auto &b:bodyForce.data()) {
        std::cout << b << std::endl;
    }
    std::cout << std::endl;
    //Parameters for compressibility
    //std::valarray<lbBase_t> alpha = LT::w0*inputAsValarray<lbBase_t>(input["fluid"]["alpha"]);
    //std::valarray<lbBase_t> alpha = LT::w0*input["fluid"]["alpha"].valarray<lbBase_t>();
    //std::valarray<lbBase_t> alpha = LT::w0*input["fluid"]["alpha"];
    std::valarray<lbBase_t> alpha = input["fluid"]["alpha"];
    alpha *= LT::w0;
    std::cout << "alpha =";
    for (int n = 0; n < nFluidFields; ++n)
        std::cout << " " << alpha[n];
    std::cout << std::endl;
    // std::valarray<lbBase_t> Gamma0 = inputAsValarray<lbBase_t>(input["fluid"]["alpha"]);
    std::valarray<lbBase_t> Gamma0 = input["fluid"]["alpha"];
    std::valarray<lbBase_t> GammaNonZero = (1-LT::w0*Gamma0)/(1-LT::w0);    
    // Phase separation: fluids
    // std::valarray<lbBase_t> beta = inputAsValarray<lbBase_t>(input["fluid"]["phaseseparation"]["beta"]);
    std::valarray<lbBase_t> beta = input["fluid"]["phaseseparation"]["beta"];
    std::cout << "beta = " << std::endl;
    for (int i=0; i < nFluidFields; ++i) {
        for (int j=0; j < nFluidFields; ++j) {
            std::cout << " " << beta[i*nFluidFields + j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    // Surface tension
    // std::valarray<lbBase_t> sigma = inputAsValarray<lbBase_t>(input["fluid"]["phaseseparation"]["sigma"]);
    std::valarray<lbBase_t> sigma = input["fluid"]["phaseseparation"]["sigma"];
    std::cout << "sigma = " << std::endl;
    for (int i=0; i < nFluidFields; ++i) {
        for (int j=0; j < nFluidFields; ++j) {
            std::cout << " " << sigma[i*nFluidFields + j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Diffusion
    // *********
    // Number of diffusive fields
    const int nDiffFields = input["diffusion"]["numfields"];
    std::cout << "nDiffFields = " << nDiffFields << std::endl;
    // Diffusion coefficients
    // std::valarray<lbBase_t> diff_coef = inputAsValarray<lbBase_t>(input["diffusion"]["fluidinteraction"]["diffcoef"]);
    std::valarray<lbBase_t> diff_coef = input["diffusion"]["fluidinteraction"]["diffcoef"];
    std::valarray<lbBase_t> diff_coef_inv(nDiffFields*nFluidFields);
    std::cout << "Diff. coef = " << std::endl;
    for (int i=0; i < nDiffFields; ++i) {
      std::cout << "diff.field[" << i <<"]: ";
        for (int j=0; j < nFluidFields; ++j) {
	    diff_coef_inv[i*nFluidFields + j] = 1./diff_coef[i*nFluidFields + j];
            std::cout << " " << diff_coef[i*nFluidFields + j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    // std::valarray<lbBase_t> betaDiff = inputAsValarray<lbBase_t>(input["diffusion"]["fluidinteraction"]["beta"]);
    std::valarray<lbBase_t> betaDiff = input["diffusion"]["fluidinteraction"]["beta"];
    std::cout << "betaDiff = " << std::endl;
    for (int k=0; k < nDiffFields; ++k) {
      std::cout << "diff.field[" << k <<"]: "<< std::endl;
      for (int i=0; i < nFluidFields; ++i) {
        for (int j=0; j < nFluidFields; ++j) {
	      std::cout << " " << betaDiff[k*nFluidFields*nFluidFields+i*nFluidFields + j];
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    // std::string dirNum = std::to_string(static_cast<int>(input["out"]["directoryNum"]));
    std::string dirNum = input["out"]["directoryNum"];
    
    std::string outputDir2 = outputDir + "/out" + dirNum;

    
    // ******************
    // MACROSCOPIC FIELDS
    // ******************
    // Fluid Flow
    // **********
    // Density
    ScalarField rho(nFluidFields, grid.size());
    ScalarField rhoRel(nFluidFields, grid.size());
    ScalarField rhoTot(1, grid.size());
    // Initiate density from file
    setScalarAttribute(rho, "init_rho_", vtklb);

    // Wall wettability
    setScalarAttributeWall(rhoRel, "init_rho_wall_", vtklb, nodes);

    // -- scale wettability
    normelizeScalarField(rhoRel, solidBoundaryNodes);

    // Velocity
    VectorField<LT> vel(nFluidFields, grid.size());
    // Initiate velocity
    for (auto fieldNo=0; fieldNo < vel.num_fields(); ++fieldNo) {
        for (auto nodeNo: bulkNodes) {
            vel.set(fieldNo, nodeNo) = 0;
        }
    }
    
    // Advection-Diffusion
    // *******************
    ScalarField phi(nDiffFields, grid.size());
    // Initiate diffusive fields
   
    for (auto nodeNo: bulkNodes) {
        for (int fieldNo=0; fieldNo < phi.num_fields(); ++fieldNo) {
    	    phi(fieldNo, nodeNo) = 0.01*rho(fieldNo, nodeNo);
	    }
        phi(2, nodeNo) = 0.01*rho(1, nodeNo);
        phi(3, nodeNo) = 0.01*rho(1, nodeNo);
        phi(4, nodeNo) = 0.01*rho(1, nodeNo);
    }
    
    // ******************
    // SETUP BOUNDARY
    // ******************
    HalfWayBounceBack<LT> bounceBackBnd(findFluidBndNodes(nodes), nodes, grid);

    // *********
    // LB FIELDS
    // *********
    // Fluid Fields
    // ************
    LbField<LT> f(nFluidFields, grid.size()); 
    LbField<LT> fTmp(nFluidFields, grid.size());
    // initiate lb distributions
    for (int fieldNo=0; fieldNo < f.num_fields(); ++fieldNo) {
        for (auto nodeNo: bulkNodes) {
            for (int q = 0; q < LT::nQ; ++q) {
                f(fieldNo, q, nodeNo) = LT::w[q]*rho(fieldNo, nodeNo);
                fTmp(fieldNo, q, nodeNo) = 0;
            }
        }
    }
    // Diffusive Fields
    //*****************
    LbField<LT> g(nDiffFields, grid.size()); 
    LbField<LT> gTmp(nDiffFields, grid.size());
    // initiate lb distributions
    for (int fieldNo=0; fieldNo < g.num_fields(); ++fieldNo) {
        for (auto nodeNo: bulkNodes) {
            for (int q = 0; q < LT::nQ; ++q) {	      
	        g(fieldNo, q, nodeNo) = LT::w[q]*phi(fieldNo, nodeNo);
                gTmp(fieldNo, q, nodeNo) = 0;
            }
        }
    }
    
    
    // **********
    // OUTPUT VTK
    // **********
    Output<LT> output(grid, bulkNodes, outputDir2, myRank, nProcs);
    output.add_file("lb_run");
    output.add_scalar_variables({"rho", "phi"}, {rho, phi});
    output.add_vector_variables({"vel"}, {vel});
    // VTK::Output<VTK_CELL, double> output(VTK::BINARY, grid.getNodePos(bulkNodes), outputDir2, myRank, nProcs);
    // output.add_file("lb_run");
    // for (int fieldNo=0; fieldNo < nFluidFields; ++fieldNo) {
    //     output.add_variable("rho" + std::to_string(fieldNo), 1, rho.get_data(), rho.get_field_index(fieldNo, bulkNodes));
    // }
    // output.add_variable("vel", LT::nD, vel.get_data(), vel.get_field_index(0, bulkNodes));
    // for (int fieldNo=0; fieldNo < nDiffFields; ++fieldNo) {
    //     output.add_variable("phi" + std::to_string(fieldNo), 1, phi.get_data(), phi.get_field_index(fieldNo, bulkNodes));
    // }
        
    // *********
    // MAIN LOOP
    // *********
    // Just ad hoc helper fields (Should be set outside the loop structure)
    std::valarray<lbBase_t> cNormInv(LT::nQ);
    std::valarray<lbBase_t> wAll(LT::nQ);
    for (int q=0; q < LT::nQNonZero_; ++q) {
        cNormInv[q] = 1.0/LT::cNorm[q];
        wAll[q] = LT::w[q];
    }
    cNormInv[LT::nQNonZero_] = 0;
    wAll[LT::nQNonZero_] = LT::w[LT::nQNonZero_];
    for (int i = 0; i <= nIterations; i++) {
        // Macroscopic values : rho, rhoTot and rhoRel = rho/rhoTot
        calcDensityFields(rho, rhoRel, rhoTot, bulkNodes, f);

        mpiBoundary.communciateScalarField(rhoRel);

        for (auto nodeNo: bulkNodes) {
            // Cacluate gradient
            ScalarField rhoRelNode(1, nFluidFields);
            for (int fieldNo = 0; fieldNo < nFluidFields; ++fieldNo)
                rhoRelNode(0, fieldNo) = rhoRel(fieldNo, nodeNo);

            const CGAttributes<LT> cgat(nFluidFields, nodeNo, cNormInv, Gamma0, rhoRelNode, rhoRel, grid);

            /* VectorField<LT> F(1, (nFluidFields*(nFluidFields-1))/2);
            ScalarField FSquare(1, F.size());
            ScalarField FNorm(1, F.size());
            LbField<LT> cDotFRC(1, F.size());
            LbField<LT> cosPhi(1, F.size());
            VectorField<LT> gradNode(1, nFluidFields);   
            int cnt = 0;
            //total effect of modified compressibility
            lbBase_t Gamma0TotNode = 0.0;
            lbBase_t GammaNonZeroTotNode = 0.0;
            
            for (int fieldNo_k=0; fieldNo_k<nFluidFields; ++fieldNo_k) {
                Gamma0TotNode += rhoRelNode(0, fieldNo_k)*Gamma0[fieldNo_k];
                GammaNonZeroTotNode += rhoRelNode(0, fieldNo_k)*GammaNonZero[fieldNo_k];
		        gradNode.set(0, fieldNo_k) = grad<LT>(rhoRel, fieldNo_k, nodeNo, grid);
                for (int fieldNo_l = 0; fieldNo_l < fieldNo_k; ++fieldNo_l) {
                    F.set(0, cnt) = rhoRelNode(0, fieldNo_l)*gradNode(0, fieldNo_k) - rhoRelNode(0, fieldNo_k)*gradNode(0, fieldNo_l);
                    cDotFRC.set(0, cnt) = LT::cDotAll(F(0,cnt));
                    FSquare(0, cnt) = LT::dot(F(0, cnt), F(0, cnt));
                    FNorm(0,cnt) = sqrt(FSquare(0, cnt));
                    if (std::abs(FNorm(0, cnt)) < lbBaseEps)
                        FNorm(0, cnt) = lbBaseEps;
                    cosPhi.set(0, cnt) = cDotFRC(0, cnt)*cNormInv/FNorm(0,cnt);
                    cnt++;                    
                }
            } */

            std::valarray<lbBase_t> feqTotRel0Node(LT::nQ);
            for (int q=0; q < LT::nQNonZero_; ++q) {
                // feqTotRel0Node[q]= wAll[q]* GammaNonZeroTotNode; //cgat.GammaNonZeroTotNode_;
                feqTotRel0Node[q]= wAll[q]* cgat.GammaNonZeroTotNode_;
            }
            // feqTotRel0Node[LT::nQNonZero_] = wAll[LT::nQNonZero_]*Gamma0TotNode; //cgat.Gamma0TotNode_;
            feqTotRel0Node[LT::nQNonZero_] = wAll[LT::nQNonZero_]*cgat.Gamma0TotNode_;

                // Calculate velocity
            // Copy of local velocity distribution
            auto velNode = calcVel<LT>(f(0, nodeNo), rhoTot(0, nodeNo));
	    lbBase_t visc_inv = rhoRelNode(0, 0)*kin_visc_inv[0];
            for (int fieldNo=1; fieldNo<nFluidFields; ++fieldNo){
                velNode += calcVel<LT>(f(fieldNo, nodeNo), rhoTot(0, nodeNo));
		visc_inv += rhoRelNode(0, fieldNo)*kin_visc_inv[fieldNo];
	        }
            velNode += 0.5*bodyForce(0, 0)/rhoTot(0, nodeNo);
	    lbBase_t tauFlNode = LT::c2Inv/visc_inv + 0.5;
            // Save density and velocity for printing
            vel.set(0, nodeNo) = velNode;
                
            // BGK-collision term
            const auto u2 = LT::dot(velNode, velNode);
            const auto cu = LT::cDotAll(velNode);
            // Calculate the Guo-force correction
            const auto uF = LT::dot(velNode, bodyForce(0, 0));
            const auto cF = LT::cDotAll(bodyForce(0, 0));

            LbField<LT> fTotNode(1,1);
            fTotNode.set(0,0) = 0;
            LbField<LT> omegaRC(1, nFluidFields);
	    
            for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo) {
	        const auto fNode = f(fieldNo, nodeNo);
                const auto rhoNode = rho(fieldNo, nodeNo);
                //const auto omegaBGK = calcOmegaBGK<LT>(fNode, tauFlNode, rhoNode, u2, cu);
		const auto feqNode = calcfeq_TEST<LT>(Gamma0[fieldNo], GammaNonZero[fieldNo], rhoNode, u2, cu);
		const auto omegaBGK = calcOmegaBGK_TEST<LT>(fNode, feqNode, tauFlNode);
                const std::valarray<lbBase_t> deltaOmegaF = rhoRel(fieldNo, nodeNo) * calcDeltaOmegaF<LT>(tauFlNode, cu, uF, cF);
                LbField<LT> deltaOmegaST(1,1);
               
                // Recoloring step
                int field_k_ind = (fieldNo*(fieldNo-1))/2;
                omegaRC.set(0, fieldNo) = 0;
                deltaOmegaST.set(0 ,0) = 0;
                for (int field_l = 0; field_l < fieldNo; ++field_l) {
                    const int F_ind = field_k_ind + field_l;
		    const int sigmaBeta_ind = fieldNo*nFluidFields + field_l;
                    deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, sigma[sigmaBeta_ind], cgat.FNorm_(0, F_ind), cgat.cDotFRC_(0, F_ind)/(cgat.FNorm_(0, F_ind)+(cgat.FNorm_(0,F_ind)<lbBaseEps)));  
                    omegaRC.set(0, fieldNo) += beta[sigmaBeta_ind]*rhoRel(field_l, nodeNo)*cgat.cosPhi_(0, F_ind);
                    /* deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, sigma[sigmaBeta_ind], FNorm(0, F_ind), cDotFRC(0, F_ind)/FNorm(0, F_ind));
                    omegaRC.set(0, fieldNo) += beta[sigmaBeta_ind]*rhoRel(field_l, nodeNo)*cosPhi(0, F_ind); */
                }
                for (int field_l = fieldNo + 1; field_l < nFluidFields; ++field_l) {
                    const int field_k_ind = (field_l*(field_l-1))/2;
                    const int F_ind =  field_k_ind + fieldNo;
		    const int sigmaBeta_ind = fieldNo*nFluidFields + field_l;
                    deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, sigma[sigmaBeta_ind], cgat.FNorm_(0, F_ind), -cgat.cDotFRC_(0, F_ind)/(cgat.FNorm_(0, F_ind)+(cgat.FNorm_(0,F_ind)<lbBaseEps)));  
                    omegaRC.set(0, fieldNo) -= beta[sigmaBeta_ind]*rhoRel(field_l, nodeNo)*cgat.cosPhi_(0,F_ind);
                    /* deltaOmegaST.set(0 ,0) += calcDeltaOmegaST<LT>(tauFlNode, sigma[sigmaBeta_ind], FNorm(0, F_ind), -cDotFRC(0, F_ind)/FNorm(0, F_ind));
                    omegaRC.set(0, fieldNo) -= beta[sigmaBeta_ind]*rhoRel(field_l, nodeNo)*cosPhi(0,F_ind); */
                }

                //omegaRC.set(0, fieldNo) *= wAll*rhoNode;
                omegaRC.set(0, fieldNo) *= rhoNode*feqTotRel0Node;

                // Calculate total lb field
                fTotNode.set(0, 0) += fNode + deltaOmegaF + omegaBGK + deltaOmegaST(0, 0);
            }
	    
	    
            // Collision and propagation
            for (int fieldNo=0; fieldNo<nFluidFields; ++fieldNo) {
                fTmp.propagateTo(fieldNo, nodeNo, rhoRel(fieldNo, nodeNo)*fTotNode(0, 0) + omegaRC(0, fieldNo), grid);
            }

	    
	    // Diffusion
	   
	    
	    
	    //lbBase_t tauDiff_aveNode = LT::c2Inv/diff_ave_inv + 0.5;
	    for (int fieldNo=0; fieldNo<nDiffFields; ++fieldNo) {
	      
	      lbBase_t diff_ave_inv = 0.0;
	      for (int fluidNo=0; fluidNo<nFluidFields; ++fluidNo){
		diff_ave_inv += rhoRel(fluidNo, nodeNo)*diff_coef_inv[fieldNo*nFluidFields+fluidNo];
	      }
	      lbBase_t tauDiff_aveNode = LT::c2Inv/diff_ave_inv + 0.5;

	      const auto gNode = g(fieldNo, nodeNo);
	      const lbBase_t phiNode = calcRho<LT>(gNode);
	      // Save density and velocity for printing
	      phi(fieldNo, nodeNo) = phiNode;

	      const auto geqNode = calcfeq_TEST<LT>(Gamma0[0], GammaNonZero[0], phiNode, u2, cu);
	      const auto omegaBGK = calcOmegaBGK_TEST<LT>(gNode, geqNode, tauDiff_aveNode);
	      g.set(fieldNo, nodeNo) = gNode + omegaBGK;
	    }

	    //------------------------------------------------------
	    //Interaction between diffusive fields and fluid fields
	    //------------------------------------------------------
	    
	    LbField<LT> omegaDI(1, nDiffFields);
	    for (int fieldNo=0; fieldNo<nDiffFields; ++fieldNo) {
	      omegaDI.set(0, fieldNo) = 0;
	    }

	    lbBase_t W, W_1, W_2, W1, W2;
	    int diffPhaseInd;
	    int solutePhaseInd;
	    LbField<LT> cosPhiTmp(1, 1);
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    //Stored interface normals are the lower triangular part of the interface normal matrix
	    //and point from phase of lower phase index toward phase of higher phase index. e.g., 0->1, 0->2, 1->2 etc.
	    //Instead of changing sign on omegaDI contribution, change sign on potential W to obtain wanted interaction result.
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 0;
	    solutePhaseInd = 0;
	    W = rhoRel(0, nodeNo) - 1;   
	    cosPhiTmp.set(0, 0) = (cgat.FNorm_(0,0)*cgat.cosPhi_(0, 0) + cgat.FNorm_(0,1)*cgat.cosPhi_(0, 1))/(cgat.FNorm_(0,0)+cgat.FNorm_(0,1));
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    /*
	    diffPhaseInd = 1;
	    solutePhaseInd = 1;
	    W = rhoRel(solutePhaseInd, nodeNo) - 0.5;
	    omegaDI.set(0, diffPhaseInd) -= betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 0]*W*cgat.cosPhi_(0, 0);
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 2]*W*cgat.cosPhi_(0, 2);
	    */
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 1;
	    solutePhaseInd = 1;
	    //surfactant 0-1-interfaces, while soluble in phase 2 
	    W1 = rhoRel(0, nodeNo) - rhoRel(1, nodeNo);
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 0]*W1*cgat.cosPhi_(0, 0);
	    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 2;
	    solutePhaseInd = 1;
	    //surfactant 0-1-interfaces
	    W1 = rhoRel(0, nodeNo); // At interface 0-1, not soluble in phase 0 (positive sign since phase 0 is lowest phase in interface) 
	    W1+= -rhoRel(1, nodeNo); // At interface 0-1, not soluble in phase 1
	    W2 = -rhoRel(2, nodeNo); // At interfaces 0-2 and 1-2, not soluble in phase 2 
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 0]*W1*cgat.cosPhi_(0, 0);
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 2]*W2*cgat.cosPhi_(0, 2);
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 2]*W2*cgat.cosPhi_(0, 1);
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 3;
	    //surfactant 1-2-interfaces
	    W1 = rhoRel(1, nodeNo) - rhoRel(2, nodeNo);
	    W2 = rhoRel(0, nodeNo);
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 1*nFluidFields + 2]*W1*cgat.cosPhi_(0, 2);
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 1]*W2*cgat.cosPhi_(0, 0);
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 2]*W2*cgat.cosPhi_(0, 1);
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 4;
	    //surfactant 0-2-interfaces
	    W1 = rhoRel(0, nodeNo); // At interface 0-2, not soluble in phase 0 (positive sign since phase 0 is lowest phase in interface) 
	    W1+= - rhoRel(2, nodeNo); // At interface 0-2, not soluble in phase 2
	    W2 = - rhoRel(1, nodeNo); // At interface 0-1, not soluble in phase 1
	    
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 2]*W1*cgat.cosPhi_(0, 1);
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 0*nFluidFields + 1]*W2*cgat.cosPhi_(0, 0);
	    omegaDI.set(0, diffPhaseInd) -= betaDiff[diffPhaseInd*nFluidFields*nFluidFields + 1*nFluidFields + 2]*W2*cgat.cosPhi_(0, 2);
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 5;
	    solutePhaseInd = 0;
	    //soluble in phase 0 (or, i.e., not soluble in phase 1 or phase 2) 
	    //W = rhoRel(0, nodeNo) - 1;
	    W_1 = -rhoRel(1, nodeNo);
	    W_2 = -rhoRel(2, nodeNo);
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W_1*cgat.cosPhi_(0, 0);
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 2]*W_2*cgat.cosPhi_(0, 1);	    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 6;
	    solutePhaseInd = 0;
	    //not soluble in phase 0 
	    W = rhoRel(0, nodeNo);
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);
	    //W_1 = -rhoRel(1, nodeNo);
	    //W_2 = -rhoRel(2, nodeNo);
	    //omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W*cgat.cosPhi_(0, 0);
	    //omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 2]*W*cgat.cosPhi_(0, 1);	    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 7;
	    solutePhaseInd = 0;
	    //W = -(rhoRel(0, nodeNo)-1);
	    //surfactant midpoint phase 0 interfaces
	    W = rhoRel(0, nodeNo) - 0.5;
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);
	    //W_1 = -rhoRel(1, nodeNo);
	    //W_2 = -rhoRel(2, nodeNo);
	    //omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W*cgat.cosPhi_(0, 0);
	    //omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 2]*W*cgat.cosPhi_(0, 1);	    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    diffPhaseInd = 8;
	    solutePhaseInd = 0;
	    //surfactant on phase 0 side of interfaces
	    W = rhoRel(0, nodeNo) - 0.75;
	    omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W*cosPhiTmp(0, 0);
	    //W_1 = -rhoRel(1, nodeNo);
	    //W_2 = -rhoRel(2, nodeNo);
	    //omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 1]*W*cgat.cosPhi_(0, 0);
	    //omegaDI.set(0, diffPhaseInd) += betaDiff[diffPhaseInd*nFluidFields*nFluidFields + solutePhaseInd*nFluidFields + 2]*W*cgat.cosPhi_(0, 1);	    
	    //-----------------------------------------------------------------------------------------------------------------------------------
	    
	    //------------------------------------------------------
	    //END Interaction between diffusive fields and fluid fields
	    //------------------------------------------------------
	    
	    
	    for (int fieldNo=0; fieldNo<nDiffFields; ++fieldNo) {
	      omegaDI.set(0, fieldNo) *= wAll*phi(fieldNo, nodeNo);
	      gTmp.propagateTo(fieldNo, nodeNo, g(fieldNo, nodeNo) + omegaDI(0, fieldNo), grid);
	    }
	    
	    
	    

        } // End nodes
        // Swap data_ from fTmp to f;
        f.swapData(fTmp);  // LBfield
	g.swapData(gTmp);  // LBfield

        // *******************
        // BOUNDARY CONDITIONS
        // *******************
        // Mpi
        mpiBoundary.communicateLbField(f, grid);
        // Half way bounce back
        bounceBackBnd.apply(f, grid);

	mpiBoundary.communicateLbField(g, grid);
        // Half way bounce back
        bounceBackBnd.apply(g, grid);

        // *************
        // WRITE TO FILE
        // *************
        if ( ((i % nItrWrite) == 0)  ) {
	  /*
	    for (auto nn: bulkNodes) {
                velIO(0, 0, nn) = vel(0, 0, nn);
                velIO(0, 1, nn) = vel(0, 1, nn);
                velIO(0, 2, nn) = 0; 

            }
            output.write("lb_run", i);
	  */
	    output.write(i);
            if (myRank==0) {
                std::cout << "PLOT AT ITERATION : " << i << std::endl;
            }
        }

    } // End iterations

    MPI_Finalize();

    return 0;
}
