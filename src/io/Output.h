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
#include <array>


//=====================================================================================
//
//                            L B   O U T P U T   U N S T R U C T U R E D 
//
//=====================================================================================
// LBOutputUnstructured is a convenience wrapper around VTK::OutputUnstructured. Template parameters:
// - T controls the output data type (double or float).
// - FMT controls ASCII vs binary (VTK::BINARY by default).
// Use LBOutputUnstructured<LT, float> to write Float32 binary for all variables.
template <typename LT, typename T=double, int FMT=VTK::BINARY, typename CELL=VTK::voxel>
class LBOutputUnstructured 
{
    private:
    VTK::OutputUnstructured<CELL,LT::nD,T> out_;
    const std::vector<int>* nodes_ = nullptr;

    public:
    //                                     LBOutputUnstructured
    //-----------------------------------------------------------------------------------
    // Construct output from a position list with an initial variable.
    LBOutputUnstructured(const std::vector<int>& pos, const std::string& dir, int rank, int nproc, const std::string& varname, const std::vector<T>& var) 
        : out_(FMT, pos, dir, rank, nproc), nodes_(nullptr)
    //-----------------------------------------------------------------------------------
    { 
        out_.add_file(varname);
        out_.add_variable(varname, var); 
    }

    //                                     LBOutputUnstructured
    //-----------------------------------------------------------------------------------
    // Construct output from a position list.
    LBOutputUnstructured(const std::vector<int>& pos, const std::string& dir, int rank, int nproc) 
        : out_(FMT, pos, dir, rank, nproc), nodes_(nullptr) { }
    //-----------------------------------------------------------------------------------

    //                                     LBOutputUnstructured
    //-----------------------------------------------------------------------------------
    // Construct output from grid nodes (common LB usage).
    LBOutputUnstructured(const Grid<LT>& grid, std::vector<int>& nodes, const std::string& dir, int rank, int nproc) 
        : out_(FMT, grid.pos(nodes), dir, rank, nproc), nodes_(&nodes) { }
    //-----------------------------------------------------------------------------------

    //                                     LBOutputUnstructured
    //-----------------------------------------------------------------------------------
    void write(double t=0.0) { out_.write(t); }
    //-----------------------------------------------------------------------------------

    //                                     LBOutputUnstructured
    //-----------------------------------------------------------------------------------
    void add_file(const std::string& name) { out_.add_file(name); }
    //-----------------------------------------------------------------------------------


    //                                     LBOutputUnstructured
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

    //                                     LBOutputUnstructured
    //-----------------------------------------------------------------------------------
    void add_vector_variables(const std::vector<std::string>& names, const std::vector<std::reference_wrapper<const VectorField<LT>>>& fields) 
    //-----------------------------------------------------------------------------------
    {
        add_variables<VectorField<LT>>(names, fields);
    }

    //                                     LBOutputUnstructured
    //-----------------------------------------------------------------------------------
    void add_scalar_variables(const std::vector<std::string>& names, const std::vector<std::reference_wrapper<const ScalarField>>& fields) 
    //-----------------------------------------------------------------------------------
    {
        add_variables<ScalarField>(names, fields);
    }

    //                                     LBOutputUnstructured
    //-----------------------------------------------------------------------------------
    template <typename F>
    void add_variable(const std::string& name, const F& field) 
    //-----------------------------------------------------------------------------------
    { 
        add_variables<F>({name}, {field});
    }
    
    //                                     LBOutputUnstructured
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

    //                                     LBOutputUnstructured
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
    //                                     LBOutputUnstructured
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

//=====================================================================================
//
//                               L B  O U T P U T  I M A G E
//
//=====================================================================================

template <typename LT, typename T=double, int FMT=VTK::BINARY>
class LBOutputImage
{
    private:
    VTK::OutputImage<LT::nD,T> out_;
    // LB-specific mapping from grid nodes to ImageData ordering.
    std::array<int, 6> piece_extent_ = {0, 0, 0, 0, 0, 0};
    std::array<int, LT::nD> origin_offset_ = {};
    std::vector<int> ordered_nodes_;

    public:
    //                                     LBOutputImage
    //-----------------------------------------------------------------------------------
    // Construct output using grid positions and validate image layout.
    LBOutputImage(const Grid<LT>& grid, const std::vector<int>& nodes, const std::string& dir, int rank, int nproc) 
        : out_(FMT, grid, nodes, dir, rank, nproc)
    //-----------------------------------------------------------------------------------
    {
        // Derive the local node bounds so ordering matches the actual block on this rank.
        std::array<int, LT::nD> min_pos;
        std::array<int, LT::nD> max_pos;
        const std::vector<int> pos = grid.pos(nodes);
        util::find_min_max<LT::nD>(pos, min_pos, max_pos);
        std::array<int, LT::nD> min_shift;
        std::array<int, LT::nD> max_shift;
        shift_bounds(min_pos, max_pos, min_shift, max_shift);
        util::build_piece_extent<LT::nD>(min_shift, max_shift, piece_extent_);
        // Precompute the node order for contiguous ImageData output.
        build_ordered_nodes(grid);
        // Sanity check that the local nodes form a full rectangular block.
        const long long expected = util::block_size<LT::nD>(min_pos, max_pos, min_shift);
        if (ordered_nodes_.size() != static_cast<size_t>(expected)) {
            std::cerr << "ERROR in LBOutputImage: Node list does not form a full rectangular block" << std::endl;
            std::cerr << "Expected " << expected << " nodes in the bounding box, got " << ordered_nodes_.size() << std::endl;
            util::safe_exit(-1);
        }
    }

    //                                     LBOutputImage
    //-----------------------------------------------------------------------------------
    void write(double t=0.0) { out_.write(t); }
    //-----------------------------------------------------------------------------------

    //                                     LBOutputImage
    //-----------------------------------------------------------------------------------
    void add_file(const std::string& name) { out_.add_file(name); }
    //-----------------------------------------------------------------------------------

    //                                     LBOutputImage
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

    //                                     LBOutputImage
    //-----------------------------------------------------------------------------------
    void add_vector_variables(const std::vector<std::string>& names, const std::vector<std::reference_wrapper<const VectorField<LT>>>& fields) 
    //-----------------------------------------------------------------------------------
    {
        add_variables<VectorField<LT>>(names, fields);
    }

    //                                     LBOutputImage
    //-----------------------------------------------------------------------------------
    void add_scalar_variables(const std::vector<std::string>& names, const std::vector<std::reference_wrapper<const ScalarField>>& fields) 
    //-----------------------------------------------------------------------------------
    {
        add_variables<ScalarField>(names, fields);
    }

    //                                     LBOutputImage
    //-----------------------------------------------------------------------------------
    template <typename F>
    void add_variable(const std::string& name, const F& field) 
    //-----------------------------------------------------------------------------------
    { 
        add_variables<F>({name}, {field});
    }
    
    //                                     LBOutputImage
    //-----------------------------------------------------------------------------------    //
    // Give a multi-field variable explicitly different names (other than just "field0", "field1", etc.) 
    // The names vector must match the number of fields.
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

    //                                     LBOutputImage
    //-----------------------------------------------------------------------------------
    void add_variables(const std::vector<std::string>& names, const std::vector<std::reference_wrapper<std::vector<T>>>& vectors) 
    //-----------------------------------------------------------------------------------
    { 
        for (auto n=names.begin(), v=vectors.begin(); n!=names.end(); ++n, ++v)
            out_.add_variable(*n, *v);
    }

    private:
    //                                     LBOutputImage
    //-----------------------------------------------------------------------------------
    template <typename F>
    void add_variable__(int i, const std::string& name, const F& field) 
    //-----------------------------------------------------------------------------------
    { 
        const int n_cells = static_cast<int>(ordered_nodes_.size());

        std::vector<int> ind(n_cells*field.dim()); 
        int n = 0;
        for (const auto node : ordered_nodes_) {
            for (auto d=0; d<field.dim(); ++d) {
                ind[n+d] = field.index(i, d, node);
            }
            n += field.dim();
        }

        auto varname { name };
        if (field.num_fields() > 1)
            varname += std::to_string(i);
        out_.add_variable(varname, field.data(), ind);
    }

    //                                     LBOutputImage
    //-----------------------------------------------------------------------------------
    // Precompute node order for the piece extent.
    void build_ordered_nodes(const Grid<LT>& grid)
    //-----------------------------------------------------------------------------------
    {
        const auto order = util::linear_order_from_extent<LT::nD>(piece_extent_);
        ordered_nodes_.clear();
        ordered_nodes_.reserve(order.size());
        for (const auto& idx : order) {
            std::array<int, LT::nD> pos;
            for (int d = 0; d < LT::nD; ++d)
                pos[d] = idx[d] + origin_offset_[d];
            ordered_nodes_.push_back(grid.nodeNo(pos));
        }
    }

    //                                     LBOutputImage
    //-----------------------------------------------------------------------------------
    // Shift local bounds into the ImageData index space and cache the origin offset.
    void shift_bounds(const std::array<int, LT::nD>& min_pos, const std::array<int, LT::nD>& max_pos, std::array<int, LT::nD>& min_shift, std::array<int, LT::nD>& max_shift)
    //-----------------------------------------------------------------------------------
    {
        const auto& origin = out_.origin();
        for (int d = 0; d < LT::nD; ++d) {
            origin_offset_[d] = static_cast<int>(origin[d]);
            min_shift[d] = min_pos[d] - origin_offset_[d];
            max_shift[d] = max_pos[d] - origin_offset_[d];
        }
    }

};

#endif /* SRC_OUTPUT_H_ */
