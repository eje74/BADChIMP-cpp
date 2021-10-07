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
// #include "../lbsolver/LBgrid.h"
// #include "../lbsolver/LBnodes.h"
// #include "../lbsolver/LBfield.h"
// #include "../lbsolver/LBvtk.h"
#include "VTK.h"


// class Output : public VTK::Output
// {
//     public:
//     Output()
//    VTK::Output<VTK::voxel, double> out_;

//    Output (std::vector<Grid>& grid) 
//     : quad_(grid[0].quad_points(), grid[0].quad_connectivity(), "out", 0, 1, VTK::BINARY), 
//       line_(grid[0].line_points(), grid[0].line_connectivity(), "out", 0, 1, VTK::BINARY), 
//       grid_(grid) 
//   {
//     //std::cout << "OUT" << std::endl;
//     int length = quad_.num_cells();
//     for (auto& g : grid)
//       g.data_vtk.resize(length*g.data_name.size());
//     quad_.add_file("streamtube");
//     line_.add_file("streamline");
//     for (auto& g : grid) {
//       int n=0;
//       for (auto& name : g.data_name) {
//         int offset = length*n++;
//         quad_.add_variable(name, g.data(), 1, offset, length);
//         line_.add_variable(name, g.data(), 1, offset, length);
//       }
//     }
//   }
//   void write(double time) {
//     for (auto& g : grid_.get()) {
//       g.update_data();
//     }
//     quad_.write(0, time);  
//     line_.write(0, time);  
//   }  

//   void calculate_points_and_connectivity(std::vector<std::vector<int>>& node_pos) {
//     // set 2D or 3D cell
//     std::vector<std::vector<int>> cell_points = (dim_.size()>2)? cell_points_3D_ : cell_points_2D_;
//     num_cell_points_ = cell_points.size();
//     connectivity_.reserve(cell_points.size()*node_pos.size());

//     std::vector<int> point_index(prod(dim_+1), -1);
//     for (const auto& n:node_pos) {
//       for (const auto& c:cell_points) {
//         std::vector<int> p = n + c;
//         // give cell-points a unique index
//         int idx = p[0] + p[1]*dim_[0] + p[2]*dim_[0]*dim_[1];
    
//         if (point_index[idx]<0) {
//           points_.insert(points_.end(), std::begin(p), std::end(p));
//           point_index[idx] = int(points_.size()/dim_.size())-1;
//         }
//         connectivity_.push_back(point_index[idx]);
//       }
//     }
//     // offset cell (corner) points by half the cell-size
//     points_ = points_ - 0.5*cell_edge_;
//   }

// };


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
    output.add_variable(fieldName, 1, val.get_data(), val.get_field_index(0, allNodes)); 
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
