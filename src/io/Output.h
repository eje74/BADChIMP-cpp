/*
 * output.h
 *
 *  Created on: 27. nov. 2017
 *      Author: janlv
 */

#ifndef SRC_OUTPUT_H_
#define SRC_OUTPUT_H_
// #include <unordered_map>
// #include <cstdio>
// #include <string>
// #include <iostream>
// #include <sstream>
// #include <fstream>
// #include <iomanip>
// #include <vector>
// #include <bitset>
// //#include <sys/types.h>
// #include <sys/stat.h>
// #if defined (_WIN32)
// #include <direct.h>
// #else
// #include <unistd.h>
// #endif
// //#include "Mpi_class.h"
// #include <mpi.h>
// #include "../lbsolver/Geo.h"
// //#include "global.h"
// #include "../lbsolver/defines.h"
// //#include "vector_func.h"
#include "../lbsolver/LBgrid.h"
#include "../lbsolver/LBnodes.h"
#include "../lbsolver/LBfield.h"
#include "../lbsolver/LBvtk.h"
#include "VTK.h"


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
    auto node_pos_all = grid.getNodePos(vtklb.beginNodeNo(), vtklb.endNodeNo());
    //auto globalDim = vtklb.getGlobaDimensions(); // Set as default to 3 dimensions, as prescribed by the Output class

    // Setup allNodes
    VTK::Output<VTK::voxel, double> output(VTK::BINARY, node_pos_all, outputDir, myRank, nProcs);

    ScalarField val(1, grid.size());
    std::vector<int> allNodes(vtklb.endNodeNo()-vtklb.beginNodeNo());
    int cnt = 0;
    for (int nodeNo = vtklb.beginNodeNo(); nodeNo < vtklb.endNodeNo(); ++nodeNo) {
        val(0, nodeNo) = scalarField[nodeNo];
        allNodes[cnt] = nodeNo;
        cnt++;
    }

    // Write field to file
    output.add_file(fieldName);
    output.add_buffer(val.get_data());
    output.add_variable(fieldName, 1, output.buffers().back().data(), val.get_field_index(0, allNodes)); 
    output.write(0);
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
