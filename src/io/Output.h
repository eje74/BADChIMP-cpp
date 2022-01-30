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
    void add_variable__(int i, const std::string& name, const Field& field, const std::vector<int>& nodes) 
    //-----------------------------------------------------------------------------------
    { 
        // Create index-vector for the field
        std::vector<int> ind(nodes.size()*field.dim()); 
        int n = 0;
        for (const auto& node : nodes) {
            for (auto d=0; d<field.dim(); ++d) {
                ind[n+d] = field.index(i, d, node);
            }
            n += field.dim();
        }
        // Add variable
        this->add_variable(name, field.dim(), field.data(), ind);
    }

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
        for (const auto& name : names ) {
            for (int i=0; i < field.num_fields(); ++i) {
                add_variable__(i, name, field, nodes);
            }
        }
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variables(const std::initializer_list<std::pair<const std::string, const Field* >>& map, const std::vector<int>& nodes) 
    //-----------------------------------------------------------------------------------
    { 
        for (const auto& [varname, field_ptr] : map) {
            int nF = field_ptr->num_fields();
            for (int i=0; i < nF; ++i) {
                // Append number to variable name if more than one field
                auto name = varname;
                if (nF > 1)
                    name += std::to_string(i);
                add_variable__(i, name, *field_ptr, nodes);
            }
        }
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_variables(const std::vector<std::string>& names, const std::vector<Field>& fields) 
    //-----------------------------------------------------------------------------------
    {
        for (size_t i=0; i<fields.size(); ++i)
            add_variables({names[i]}, fields[i]);    
    } 
};

#endif /* SRC_OUTPUT_H_ */
