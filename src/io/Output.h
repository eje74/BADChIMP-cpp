/*
 * output.h
 *
 *  Created on: 27. nov. 2017
 *      Author: janlv
 */

#ifndef SRC_OUTPUT_H_
#define SRC_OUTPUT_H_
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBgrid.h"
#include "../lbsolver/LBfield.h"
#include "VTK.h"


// Trick to choose pixel or voxel depending on the dimension
template <bool B>
struct CELL { typedef VTK::voxel type; };
template <>
struct CELL<false> { typedef VTK::pixel type; };


//=====================================================================================
//
//                                    O U T P U T 
//
//=====================================================================================

template <typename LT, typename T=double, typename C=typename CELL<(LT::nD>2)>::type>
class Output : public VTK::Output<C,T>
{
    public:
    //                                     Output
    //-----------------------------------------------------------------------------------
    Output(const std::vector<int>& pos, const std::string& dir, int rank, int nproc, const std::string& varname, const std::vector<T>& var, int dim=1) 
        : VTK::Output<C,T>(VTK::BINARY, pos, dir, rank, nproc) 
    //-----------------------------------------------------------------------------------
    { 
        this->add_file(varname);
        this->add_variable(varname, dim, var); 
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    Output(const std::vector<int>& pos, const std::string& dir, int rank, int nproc) 
        : VTK::Output<C,T>(VTK::BINARY, pos, dir, rank, nproc) { }
    //-----------------------------------------------------------------------------------


    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variables(const std::vector<std::string>& names, const Field& field, const std::vector<int>& nodes) 
    //-----------------------------------------------------------------------------------
    { 
        auto a = names.size();
        auto b = field.num_fields();
        if (int(a) != b) {
            std::cerr << "*** ERROR in Output::add_variables: " << a << " variable names are given for " << b << " fields" << std::endl;
            exit(1); 
        }
        std::vector<int> ind(nodes.size()*field.dim()); 
        for (int i=0; i < field.num_fields(); ++i) {
            int n = 0;
            for (const auto& node : nodes) {
                for (auto d=0; d<field.dim(); ++d) {
                    ind[n+d] = field.index(i, d, node);
                }
                n += field.dim();
            }
            this->add_variable(names[i], field.dim(), field.data(), ind);
        }
    }

};



// /* ********************************************************************** *
//  *                                                                        *
//  *          FUNCTIONS                                                     *
//  *                                                                        *
//  *                                                                        *
//  * ********************************************************************** */

// template<typename DXQY, typename T>
// void outputStdVector(const std::string &fieldName, const std::vector<T> &scalarField, const std::string &outputDir, 
//     const int &myRank, const int &nProcs, const Grid<DXQY> &grid, const LBvtk<DXQY> &vtklb)
// {

//     // Write field to file
//     VTK::Output<typename CELL<(DXQY::nD>2)>::type, int> output(VTK::BINARY, grid.pos(), outputDir, myRank, nProcs);
//     output.add_file(fieldName);
//     output.add_variable(fieldName, 1, scalarField); 
//     output.write(0);
// }

// template<typename DXQY>
// void outputGeometry(const std::string &fileName, const std::string &outputDir, const int &myRank, const int &nProcs, 
//     const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, const LBvtk<DXQY> &vtklb)
// {
    
//     std::vector<int> val(grid.size());
//     val[0] = -1;
//     for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
//         val[nodeNo] = nodes.isSolid(nodeNo) ? 1 : 0;        
//     }
    
//     outputStdVector(fileName, val, outputDir, myRank, nProcs, grid, vtklb);

// }


#endif /* SRC_OUTPUT_H_ */
