/*
 * output.h
 *
 *  Created on: 27. nov. 2017
 *      Author: janlv
 */

#ifndef SRC_OUTPUT_H_
#define SRC_OUTPUT_H_
//#include <type_traits>
//#include "../lbsolver/LBgrid.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBgrid.h"
#include "../lbsolver/LBfield.h"
//#include "../lbsolver/LBvtk.h"
#include "VTK.h"

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
    private:
    const std::vector<int>& bulk_nodes_;

    public:
    //                                     Output
    //-----------------------------------------------------------------------------------
    Output(const Grid<LT>& grid, const std::vector<int>& bulk_nodes, const std::string& dir, int rank, int nproc) 
        : VTK::Output<C,T>(VTK::BINARY, grid.getNodePos(bulk_nodes), dir, rank, nproc), bulk_nodes_(bulk_nodes) { }
    //-----------------------------------------------------------------------------------

    // //                                     Output
    // //-----------------------------------------------------------------------------------
    // Output(const std::vector<int>& node_pos, const std::string& dir, int rank, int nproc) 
    //     : VTK::Output<C,T>(VTK::BINARY, node_pos, dir, rank, nproc), bulk_nodes_(bulk_nodes) { }
    // //-----------------------------------------------------------------------------------

    //                                     Output
    //-----------------------------------------------------------------------------------
    void check_size(int a, int b)
    //-----------------------------------------------------------------------------------
    {
        if (a != b) {
            std::cerr << "*** ERROR in add_variables: names.size() (" << a << ") != field.size() (" << b << ")" << std::endl;
            exit(1); 
        }
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    std::vector<int> index(int field, const ScalarField& scalar) const 
    //-----------------------------------------------------------------------------------
    {
        std::vector<int> ind;
        ind.reserve(bulk_nodes_.size());
        const auto* begin = &scalar.data()[0];
        for (const auto& node : bulk_nodes_) {
            ind.push_back( &scalar(field, node) - begin );
        }
        return ind;
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    std::vector<int> index(int field, const VectorField<LT>& vector) const 
    //-----------------------------------------------------------------------------------
    {
      std::vector<int> ind; 
      ind.reserve(bulk_nodes_.size()*LT::nD);
      const auto* begin = &vector.data()[0];
      for (const auto& node : bulk_nodes_) {
        for (auto d=0; d<LT::nD; ++d) {
          ind.push_back( &vector(field, d, node) - begin );
        }
      }
      return ind;
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variables(const std::vector<std::string>& names, const ScalarField& sfield) 
    //-----------------------------------------------------------------------------------
    { 
        check_size(names.size(), sfield.num_fields());
        for (int i=0; i < sfield.num_fields(); ++i)
            VTK::Output<C,T>::add_variable(names[i], 1, sfield.data(), index(i, sfield));
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variables(const std::vector<std::string>& names, const VectorField<LT>& vfield) 
    //-----------------------------------------------------------------------------------
    { 
        check_size(names.size(), vfield.num_fields());
        for (int i=0; i < vfield.num_fields(); ++i)
            VTK::Output<C,T>::add_variable(names[i], LT::nD, vfield.data(), index(i, vfield));
    }
};


/* ********************************************************************** *
 *                                                                        *
 *          FUNCTIONS                                                     *
 *                                                                        *
 *                                                                        *
 * ********************************************************************** */

template<typename DXQY, typename T>
void outputStdVector(const std::string &fieldName, const std::vector<T> &scalarField, const std::string &outputDir, 
    const int &myRank, const int &nProcs, const Grid<DXQY> &grid, const LBvtk<DXQY> &vtklb)
{

    // ScalarField val(1, grid.size());
    // std::vector<int> allNodes(vtklb.endNodeNo()-vtklb.beginNodeNo());
    // int cnt = 0;
    // for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
    //     val(0, nodeNo) = scalarField[nodeNo];
    //     allNodes[cnt] = nodeNo;
    //     std::cout << nodeNo << ", "<< std::endl;
    //     cnt++;
    // }

    // Write field to file
    VTK::Output<typename CELL<(DXQY::nD>2)>::type, int> output(VTK::BINARY, grid.all_nodes(), outputDir, myRank, nProcs);
    //Output<DXQY> output(grid, grid.all_nodes(), outputDir, myRank, nProcs);
    output.add_file(fieldName);
    //output.add_variables({fieldName}, val); 
    output.add_variable(fieldName, 1, scalarField); 
    output.write(0);
    // VTK::Output<VTK::voxel, double> output(VTK::BINARY, grid.getNodePos(vtklb.beginNodeNo(), vtklb.endNodeNo()), outputDir, myRank, nProcs);
    // output.add_file(fieldName);
    // output.add_variable(fieldName, 1, val.data(), val.index(0, allNodes)); 
    // output.write(0);
}

template<typename DXQY>
void outputGeometry(const std::string &fileName, const std::string &outputDir, const int &myRank, const int &nProcs, 
    const Nodes<DXQY> &nodes, const Grid<DXQY> &grid, const LBvtk<DXQY> &vtklb)
{
    
    std::vector<int> val(grid.size());
    val[0] = -1;
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        val[nodeNo] = nodes.isSolid(nodeNo) ? 1 : 0;        
    }
    
    outputStdVector(fileName, val, outputDir, myRank, nProcs, grid, vtklb);

}


#endif /* SRC_OUTPUT_H_ */
