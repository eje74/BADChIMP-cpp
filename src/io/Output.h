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


//=====================================================================================
//
//                                    O U T P U T 
//
//=====================================================================================

template <typename LT, typename T=double, int FMT=VTK::BINARY, typename CELL=VTK::voxel>
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
  void write(double t=0.0, bool overwrite=false) { out_.write(t, overwrite); }
    //-----------------------------------------------------------------------------------

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_file(const std::string& name) { out_.add_file(name); }
    //-----------------------------------------------------------------------------------


    //                                     Output
    //-----------------------------------------------------------------------------------
    // Variable names are given by a vector of the same size as the vector of variables.
    // For multi-field variables the name is appended by the field number.
    //
    // Examples: 
    //          add_variables({"rho","eff_nu"},{rho, eff_nu})    // Separate functions for ScalarFields and
    //          add_variables({"vel","force"},{vel, force_tot})  // VectorFields
    //-----------------------------------------------------------------------------------
    template <typename F>
    void add_variables(const std::vector<std::string>& names, const std::vector<std::reference_wrapper<const F>>& fields) 
    //-----------------------------------------------------------------------------------
    { 
        for (size_t i=0; i<names.size(); ++i) {
            for (int f=0; f < fields[i].get().num_fields(); ++f) {
                add_variable__(f, names[i], fields[i].get());
            }
        }
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_vector_variables(const std::vector<std::string>& names, const std::vector<std::reference_wrapper<const VectorField<LT>>>& fields) 
    //-----------------------------------------------------------------------------------
    {
        add_variables<VectorField<LT>>(names, fields);
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    void add_scalar_variables(const std::vector<std::string>& names, const std::vector<std::reference_wrapper<const ScalarField>>& fields) 
    //-----------------------------------------------------------------------------------
    {
        add_variables<ScalarField>(names, fields);
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    template <typename F>
    void add_variable(const std::string& name, const F& field) 
    //-----------------------------------------------------------------------------------
    { 
        add_variables<F>({name}, {field});
    }
    
    //                                     Output
    //-----------------------------------------------------------------------------------    //
    // Give a multi-field variable explicitly different names (other than just "field0", "field1", etc.) 
    // The names vector must match the number of fields.
    //
    // Example: 
    //            add_variable_names({"rho_one","rho_two"}, rho)
    //-----------------------------------------------------------------------------------
    template <typename F>
    void add_variable_with_names(const std::vector<std::string>& names, const F& field) 
    //-----------------------------------------------------------------------------------
    { 
        if (static_cast<int>(names.size()) != field.num_fields()) {
            std::cerr << "ERROR in add_variables: Size of names-vector does not match number of fields, " << names.size() << " != " << field.num_fields() << "\n";
            util::safe_exit(-1);
        }
        int i = 0;
        for (const auto& name : names) { 
            add_variable__(i++, name, field);
        }
    }

    //                                     Output
    //-----------------------------------------------------------------------------------
    // Add std::vector, not ScalarField or VectorField. Variable names are given by a vector 
    // of the same size as the vector of variable-references.
    //
    // Examples: 
    //          add_variables({"geo"}, {geo})
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
    template <typename F>
    void add_variable__(int i, const std::string& name, const F& field) 
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
        // Append number to variable name if more than one field
        auto varname { name };
        if (field.num_fields() > 1)
            varname += std::to_string(i);
        // Add variable        
        out_.add_variable(varname, field.data(), ind);
    }


};

#endif /* SRC_OUTPUT_H_ */
