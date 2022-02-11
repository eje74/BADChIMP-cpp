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

  struct lattice_3D {
    static constexpr char name[] = "D3Q19";
    static constexpr int type = VTK::vertex::type;
    static constexpr bool center = false;
    static constexpr int n = 19;
    static constexpr std::array<std::array<int,3>, n> points = { {{1,0,0}, 
                                                                  {0,1,0}, 
                                                                  {0,0,1}, 
                                                                  {1,1,0}, 
                                                                  {1,-1,0}, 
                                                                  {1,0,1}, 
                                                                  {1,0,-1}, 
                                                                  {0,1,1}, 
                                                                  {0,1,-1}, 
                                                                  {-1,0,0}, 
                                                                  {0,-1,0}, 
                                                                  {0,0,-1}, 
                                                                  {-1,-1,0}, 
                                                                  {-1,1,0}, 
                                                                  {-1,0,-1}, 
                                                                  {-1,0,1}, 
                                                                  {0,-1,-1}, 
                                                                  {0,-1,1}, 
                                                                  {0,0,0}}};
  };

//=====================================================================================
//
//                                    O U T P U T 
//
//=====================================================================================

template <typename LT, typename T=double, int FMT=VTK::ASCII, typename CELL=lattice_3D>
class Output 
{
    private:
    VTK::Output<CELL,LT::nD,T> out_;
    const std::vector<int>* nodes_ = nullptr;

    public:
    //                                     Output
    //-----------------------------------------------------------------------------------
    Output(const std::vector<int>& pos, const std::string& dir, int rank, int nproc, const std::string& varname, const std::vector<T>& var) 
        : out_(FMT, pos, dir, rank, nproc), nodes_(nullptr)
    //-----------------------------------------------------------------------------------
    { 
        out_.add_file(varname);
        out_.add_variable(varname, var); 
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    Output(const std::vector<int>& pos, const std::string& dir, int rank, int nproc) 
        : out_(FMT, pos, dir, rank, nproc), nodes_(nullptr) { }
    //-----------------------------------------------------------------------------------

    //                                     Output
    //-----------------------------------------------------------------------------------
    Output(const Grid<LT>& grid, std::vector<int>& nodes, const std::string& dir, int rank, int nproc) 
        : out_(FMT, grid.pos(nodes), dir, rank, nproc), nodes_(&nodes) { }
    //-----------------------------------------------------------------------------------

    //                                     Output
    //-----------------------------------------------------------------------------------
    void write(double t=0.0) {out_.write(t);}
    //-----------------------------------------------------------------------------------

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_file(const std::string& name) {out_.add_file(name);}
    //-----------------------------------------------------------------------------------


    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variables(const std::vector<std::string>& names, const std::vector<std::reference_wrapper<const Field>>& fields) 
    //-----------------------------------------------------------------------------------
    { 
        for (size_t i=0; i<names.size(); ++i) {
            int nF = fields[i].get().num_fields();
            for (int f=0; f < nF; ++f) {
                // Append number to variable name if more than one field
                auto name = names[i];
                if (nF > 1)
                    name += std::to_string(f);
                add_variable__(f, name, fields[i]);
            }
        }
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variable(const std::string& name, const Field& field) 
    //-----------------------------------------------------------------------------------
    { 
        add_variables({name}, {field});
    }
    
    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variables(const std::initializer_list<std::pair<std::vector<std::string>, std::reference_wrapper<Field>>>& map) 
    //-----------------------------------------------------------------------------------
    { 
        for (const auto& [names, field] : map) {
            int i = 0;
            for (const auto& name : names) { 
                add_variable__(i++, name, field);
            }
        }
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variables(const std::vector<std::string>& names, const std::vector<std::reference_wrapper<std::vector<T>>>& vectors) 
    //-----------------------------------------------------------------------------------
    { 
        for (auto n=names.begin(), v=vectors.begin(); n!=names.end(); ++n, ++v)
            out_.add_variable(*n, *v);
    }
    

    private:
    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variable__(int i, const std::string& name, const Field& field) 
    //-----------------------------------------------------------------------------------
    { 
        // Create index-vector for the field
        const auto& nodes = *nodes_;
        std::vector<int> ind(nodes.size()*field.dim()); 
        int n = 0;
        for (const auto& node : nodes) {
            for (auto d=0; d<field.dim(); ++d) {
                ind[n+d] = field.index(i, d, node);
            }
            n += field.dim();
        }
        // Add variable        
        out_.add_variable(name, field.data(), ind);
    }


};

#endif /* SRC_OUTPUT_H_ */
