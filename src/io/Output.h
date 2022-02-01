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
    private:
    const std::vector<int>* nodes_ = nullptr;

    public:
    //                                     Output
    //-----------------------------------------------------------------------------------
    Output(const std::vector<int>& pos, const std::string& dir, int rank, int nproc, const std::string& varname, const std::vector<T>& var, int dim=1) 
        : VTK::Output<C,T>(VTK::BINARY, pos, dir, rank, nproc), nodes_(nullptr)
    //-----------------------------------------------------------------------------------
    { 
        this->add_file(varname);
        this->add_variable(varname, var); 
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    Output(const std::vector<int>& pos, const std::string& dir, int rank, int nproc) 
        : VTK::Output<C,T>(VTK::BINARY, pos, dir, rank, nproc), nodes_(nullptr) { }
    //-----------------------------------------------------------------------------------

    //                                     Output
    //-----------------------------------------------------------------------------------
    Output(const Grid<LT>& grid, std::vector<int>& nodes, const std::string& dir, int rank, int nproc) 
        : VTK::Output<C,T>(VTK::BINARY, grid.pos(nodes), dir, rank, nproc), nodes_(&nodes) { }
    //-----------------------------------------------------------------------------------


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
        //this->add_variable(name, field.data(), ind);
        VTK::Output<C,T>::add_variable(name, field.data(), ind);
    }

    // //                                     Output
    // //-----------------------------------------------------------------------------------
    // //  
    // //  
    // void add_variables(const std::vector<std::string>& names, const Field& field) 
    // //-----------------------------------------------------------------------------------
    // { 
    //     auto a = names.size();
    //     auto b = field.num_fields();
    //     if (int(a) != b) {
    //         std::cerr << "*** ERROR in Output::add_variables: " << a << " variable names are given for " << b << " fields" << std::endl;
    //         exit(1); 
    //     }
    //     for (const auto& name : names ) {
    //         for (int i=0; i < field.num_fields(); ++i) {
    //             add_variable__(i, name, field);
    //         }
    //     }
    // }

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
    void add_variables(const std::initializer_list<std::pair<const std::string, const std::reference_wrapper<Field>>>& map) 
    //-----------------------------------------------------------------------------------
    { 
        for (const auto& [name, field] : map) {
            add_variables({name}, {field});
        }
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variables(const std::vector<std::string>& names, const std::vector<std::reference_wrapper<std::vector<T>>>& vectors) 
    //-----------------------------------------------------------------------------------
    { 
        for (auto n=names.begin(), v=vectors.begin(); n!=names.end(); ++n, ++v)
            this->add_variable(*n, *v);
    }
    
};

#endif /* SRC_OUTPUT_H_ */
