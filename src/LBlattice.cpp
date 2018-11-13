#include "LBlattice.h"

Lattice::Lattice(const int nDirections, const int nDimensions) : nDirections_(nDirections), nDimensions_(nDimensions)  // Lattice constructor
{
    w_ = new double [nDirections_];                     /* List of lattice weights */
    cDMajor_ = new int [nDirections_ * nDimensions_];   /*  Basis velocities on the form c[alpha][dim]
                                                         *  So that c = [c(alpha=0)_x|c(alpha=0)_y, ...]
                                                         */
    cQMajor_ = new int [nDirections_ * nDimensions_];   /* Basis velocities on the form c[dim][alpha].
                                                         * So that c = [c(alpha=0)_x|c(alpha=1)_x, ...]
                                                         */
    reverseStep_ = nDirections_ / 2;                    /* Number of bounce back pairs */
    nNonZeroDirections_ = 2 * reverseStep_;             /* Number of non-rest particles velocities */
}

Lattice::~Lattice()  // Lattice destructor
{
    delete [] w_;
    delete [] cDMajor_;
    delete [] cQMajor_;
//        std::cout << "Destructed a Lattice object" << std::endl;
}


