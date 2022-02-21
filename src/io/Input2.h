/*
 * Input.h
 *
 *  Created on: 4. mar. 2015
 *      Author: janlv
 */

#ifndef SRC_INPUT_H_
#define SRC_INPUT_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
//#include <valarray>
#include <map>
#include <algorithm>
#include <numeric>
#include <execution>


//------------------------------------------------------------
// Demonstration of the Input-class
// 
// The contents of the example.file below can be access in a C++
// program like this:
//
// #include <Input.h>
//
// Input input("example.file");
// double version = input["version"];                       // version = 1.1
// double pi = input["math"]["pi"];                         // pi = 3.12
// std::vector<int> fibo = input["math"]["series"]["fibo"]; // fibo = 2, 3, 5, 7, 11, 13, 17, 19
// std::vector<int> binary = input["binary"];               // binary = 0, 1, 0, 1, 1, 1, 1, 0
//
// 
// The contents of the example.file:
// # Use # to comment whole lines or part of a line
// #
// version 1.1                 # Both numbers and 
// dir myoutput/out1 save_dir  # text are valid input
// day_sec 1/86400             # Perform simple math operations
// time $day_sec*100           # Use day_sec as a variable
// add 2+2
// sub 1-10
// mult 2*4
// <math>                      # Make a new section using specified keyword tags
//    pi 3.14
//    e  2.72
//    <series>
//       fibo 1 1 2 3 5 8
//       <prime int>           # Unless 'int' or 'char' is specified, a 'double' datatype is assumed
//          2 3 5 7
//          11 13 17 19
//       <end>
//    <end>
// <end>
// <binary char>               # use a char datatype to read digit by digit
//    0101
//    1110
// <end>
//
//
// Numerical input-values can be used in math operations
// double sum = input["e"]+input["pi"];
// 
// String-values can be appended to other strings
// std::string str = "output/" + input["outdir"];
// outdir += input["outdir"];



namespace str_func
{

    //-----------------------------------------------------------------------------------
    // Check if a string containing multiple words can be converted to numbers
    //-----------------------------------------------------------------------------------
    bool is_numeric(const std::string& str)
    {
        std::istringstream ss(str);
        double dbl;
        while (not ss.eof()) {
            ss >> dbl;   // Try to read the number
            if (ss.fail()) { 
                return false;
            } 
        }
        return true; 
    }

}


//=====================================================================================
//
//                                    B L O C K 
//
// A block is the content between an opening and a closing tag
// A block can contain blocks recursively
// A block is accessed by the []-operator
//=====================================================================================

class Block
{
public:
    int level_ = 0;                     // Used for indenting output
    std::string name_;                  // name of the block given by the first string on the line. If no string, the line number is used
    std::vector<Block> blocks_;         // a block is recursive and can contain other (sub-) blocks
    std::vector<double> values_;        // Block num values 
    std::vector<std::string> strings_;  // Block string values
    std::string datatype_;              // Default type is double, but int and char can be specified
    Block* parent_ = nullptr;
    Block* not_found_ = nullptr;
    bool missing_ok_ = false;
    // std::array<size_t,2> pos_ = {0,0};
    std::string filename_;

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block() : name_(), blocks_(), values_(), strings_(), datatype_() { }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block(const std::string &name, Block& not_found, bool missing_ok=false) 
        : name_(name), blocks_(), values_(), strings_(), datatype_(), not_found_(&not_found), missing_ok_(missing_ok), filename_(name) { }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block(const std::string &name, int level, Block& parent) 
        : level_(level), name_(name), blocks_(), values_(), strings_(), datatype_(), parent_(&parent), 
          not_found_(parent.not_found_), missing_ok_(parent.missing_ok_), filename_(parent.filename_)
    //-----------------------------------------------------------------------------------
    {
        if ( str_func::is_numeric(name_) ) {
            // Use line-number as name. A new Block is created for each line
            name_ = "0";
            if (parent_)
                name_ = std::to_string(parent_->blocks_.size());
        }
    }

    // //                                     Block
    // //-----------------------------------------------------------------------------------
    // Block(const Block& b) : name_(b.name_), blocks_(b.blocks_), values_(b.values_), strings_(b.strings_), datatype_(b.datatype_), 
    //                         parent_(b.parent_), not_found_(b.not_found_), missing_ok_(b.missing_ok_), pos_(b.pos_), filename_(b.name) { }
    // //-----------------------------------------------------------------------------------
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    void info() const
    //-----------------------------------------------------------------------------------
    {
        std::cout << name_ << ", level = " << level_ << ", blocks = ";
        for (const auto& b : blocks_)
            std::cout << b.name_ << ", ";
        std::cout << "values = " << values_.size();
        std::cout << ", strings = " << strings_.size() << ", datatype = " << datatype_ << ", parent = " << ((parent_)? ("None") : parent_->name_);
        std::cout << std::endl;                
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block* parent() const { if (parent_) { return (parent_); } else { error("Block " + name_ + " has no parent!"); return nullptr;} }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    bool is_vector() const { return (values_.size() > 1) ? true : false ; }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    operator int() const { return static_cast<int>(values_[0]); }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    operator bool() const { return (name_.empty()) ? false : true; }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    operator double() const { return values_[0]; }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    operator std::string() const 
    //-----------------------------------------------------------------------------------
    { 
        if (strings_.empty())
            return std::to_string(values_[0]);
        else 
            return strings_[0]; 
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Arithmetic block,block operators
    //
    friend double operator+(const Block& lhs, const Block& rhs) { return lhs.values_[0]+rhs.values_[0]; }
    friend double operator-(const Block& lhs, const Block& rhs) { return lhs.values_[0]-rhs.values_[0]; }
    friend double operator*(const Block& lhs, const Block& rhs) { return lhs.values_[0]*rhs.values_[0]; }
    friend double operator/(const Block& lhs, const Block& rhs) { return lhs.values_[0]/rhs.values_[0]; }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    // String operators
    //
    friend std::string operator+(const std::string& lhs, const Block& block) { return lhs + block.strings_[0]; }
    friend std::string operator+(const Block& block, const std::string& rhs) { return block.strings_[0] + rhs; }
    friend std::string& operator+=(std::string& lhs, const Block& block) { lhs = lhs + block.strings_[0]; return lhs; }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Arithmetic number,block and block,number operators
    //
    template <typename T>
    friend double operator+(const T& lhs, const Block& block) { return block.is_vector() ? block.error("vector") : static_cast<double>(lhs) + block.values_[0]; }
    template <typename T>
    friend double operator+(const Block& block, const T& rhs) { return block.is_vector() ? block.error("vector") : block.values_[0] + static_cast<double>(rhs); }
    template <typename T>
    friend double operator-(const T& lhs, const Block& block) { return block.is_vector() ? block.error("vector") : static_cast<double>(lhs) - block.values_[0]; }
    template <typename T>
    friend double operator-(const Block& block, const T& rhs) { return block.is_vector() ? block.error("vector") : block.values_[0] - static_cast<double>(rhs); }
    template <typename T>
    friend double operator*(const T& lhs, const Block& block) { return block.is_vector() ? block.error("vector") : static_cast<double>(lhs) * block.values_[0]; }
    template <typename T>
    friend double operator*(const Block& block, const T& rhs) { return block.is_vector() ? block.error("vector") : block.values_[0] * static_cast<double>(rhs); }
    template <typename T>
    friend double operator/(const T& lhs, const Block& block) 
    { 
        if (block.is_vector())
            return block.error("vector");
        else if (block.values_[0] != 0.0) 
            return static_cast<double>(lhs)/block.values_[0];
        else
            return nanf(""); 
    }
    template <typename T>
    friend double operator/(const Block& block, const T& rhs) 
    { 
        if (block.is_vector())
            return block.error("vector");
        else if (rhs != 0) 
            return block.values_[0]/static_cast<double>(rhs);
        else
            return nanf(""); 
    }
    template <typename T>
    friend int operator%(const T& lhs, const Block& block) { return block.is_vector() ? block.error("vector") : lhs % static_cast<int>(block.values_[0]); }
    template <typename T>
    friend int operator%(const Block& block, const T& rhs) { return block.is_vector() ? block.error("vector") : static_cast<int>(block.values_[0]) % rhs; }
    //-----------------------------------------------------------------------------------
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    // Return number of rows
    size_t nrows(void) const { return blocks_.size(); }
    //-----------------------------------------------------------------------------------
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    // Return number of columns
    size_t ncols(void) const { return (values_.size()>0) ? values_.size() : strings_.size(); }
    //-----------------------------------------------------------------------------------
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    // return number of rows matching given string
    int nrows_match(const std::string &pattern) const
    //-----------------------------------------------------------------------------------
    {
        int m=0;
        for (const auto& b:blocks_) {
            // if (b->name_ == pattern)
            if (b.name_ == pattern)
            ++m;
        }
        return m;
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    int nrows_not_match(const std::string &pattern) const { return nrows()-nrows_match(pattern); }
    //-----------------------------------------------------------------------------------
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    // Implicit conversion returning a std::vector of typename T
    // Example: std::vector<int> size = input["size-vector"];
    //
    template <typename T> 
    operator std::vector<T>() const       
    //-----------------------------------------------------------------------------------
    {
        // std::cout << "operator std::vector<T>()" << ", " << name_ << std::endl;
        if (nrows()>0) {
            // If the block contains several unnamed lines, return a flattened
            // vector of all lines of the block.
            // See <prime> or <binary> in the example.file on top.
            std::vector<T> vec;
            vec.reserve(nrows()*blocks_[0].ncols());
            for (const auto& bl:blocks_) {
                vec.insert(vec.end(), bl.values_.begin(), bl.values_.end());
            }
            return vec;
        } else {
            // only one line, create typename vector and return it
            return(std::vector<T>(values_.begin(), values_.end()));
        }
    }   

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Implicit conversion returning a std::vector of std::string
    //
    operator std::vector<std::string>() const { return(std::vector<std::string>(strings_.begin(), strings_.end())); }   
    //-----------------------------------------------------------------------------------


    //                                     Block
    //-----------------------------------------------------------------------------------
    // Return the values vector
    //
    template <typename T>
    std::vector<T> row() const { return(std::vector<T>(values_.begin(), values_.end())); }
    //-----------------------------------------------------------------------------------
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    // Create and return a 1D VECTOR of column n of a block
    //
    template <typename T>
    std::vector<T> column(const int n) const   
    //-----------------------------------------------------------------------------------
    {
        std::vector<T> col;
        for (const auto& bl:blocks_) {
            if ( (n+1)>int(bl.values_.size())) {
                std::cerr << "ERROR in Block::get_column: index " << n
                    << " beyond limit " << bl.values_.size()-1 << std::endl;
                exit(1);
            }
            col.push_back(bl.values_[n]);
        }
        return col;
    }

    
    //                                     Block
    //-----------------------------------------------------------------------------------
    // Create and return a MATRIX (2D vector) of size nrows by ncols
    //
    template <typename T>
    std::vector< std::vector<T> > matrix() const           
    //-----------------------------------------------------------------------------------
    {
        std::vector< std::vector<T> > mat(nrows(), std::vector<T>(blocks_[0].ncols()));
        int n = 0;
        for (auto &row : mat) {
            row = blocks_[n++].row<T>();
        }
        return mat;
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Create and return a SYMMETRIC MATRIX of size nrows by ncols
    // where the lower triangle is overwritten by the upper triangle
    // 
    template <typename T>
    std::vector< std::vector<T> > symmetric_matrix() const
    //-----------------------------------------------------------------------------------
    {
        std::vector< std::vector<T> > mat = matrix<T>();
        int rows = mat.size();
        int cols = mat[0].size();
        if ( rows != cols) {
            std::cerr << "ERROR! The matrix is not quadratic (" << rows << ", " << cols << "), unable to create symmetric matrix!" << std::endl;
            exit(1);
        }
        for (int i=0; i<rows; ++i) {
            for (int j=0; j<cols; ++j) {
                if (j>i) {
                    mat[j][i] = mat[i][j];
                }
            }
        }
        return mat;
    }
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    std::vector<std::string> names(void) const      
    //-----------------------------------------------------------------------------------
    {
        std::vector<std::string> var_names;
        for (size_t i=0; i<blocks_.size(); ++i) {
            var_names.emplace_back(blocks_[i].name_);
        }
        return var_names;
    }
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    const Block& operator[](const char* keyword) const
    //-----------------------------------------------------------------------------------
    {
        if (blocks_.empty()) {
            return *this;
        } else {
            if ( const Block* block = find(std::string(keyword)) ) {
                return *block;
            } else {
                if (not missing_ok_) 
                    error("keyword " + std::string(keyword) + " not found!");
                return *not_found_;
            }
        }
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    double operator[](const int ind) const
    //-----------------------------------------------------------------------------------
    {
        if (int(values_.size())<ind+1) {
            std::cerr << std::endl << "ERROR reading input-file: '" << name_ << "': index outside limit, " << ind << " > " << values_.size()-1 << std::endl << std::endl;
            std::cerr << values_[0] << std::endl;
            exit(1);
        }
        return values_[ind];
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    //
    template <typename T>
    void push_back(const std::string& word) 
    //-----------------------------------------------------------------------------------
    {
        T w;
        std::istringstream iss(word);
        while (iss >> w) {
            if (typeid(w) == typeid(char))
                w -= '0'; // Get value from ascii 
            values_.push_back(w);
        }
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    void push_back_number(const std::string& word)
    //-----------------------------------------------------------------------------------
    {
        // Numeric value, add it
        auto dtype = parent()->datatype_;
        if (dtype == "char") {
            push_back<char>(word);
        } else if (dtype == "int") {
            push_back<int>(word);
        } else {
            push_back<double>(word);
        }
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    friend std::ostream& operator<<(std::ostream& out, const Block& block)
    //-----------------------------------------------------------------------------------
    {
        std::string indent = std::string(block.level_*3, ' ');
        // out << indent+block.name_ << " (" << block.pos_[0] << "," << block.pos_[1] << ")" << ": ";
        out << indent+block.name_ << ": ";
        for (const auto& val : block.values_) {
            out << val << " ";
        }
        for (const auto& str : block.strings_) {
            out << str << " ";
        }
        out << std::endl;
        // recursive call
        for (const auto& block : block.blocks_) {
            out << block;
        }
        return out;
    }


private:

    //                                     Block
    //-----------------------------------------------------------------------------------
    //
    int error(const std::string& msg) const
    //-----------------------------------------------------------------------------------
    {
        std::cerr << std::endl << "*** ERROR during reading of input-file ***" << std::endl;
        std::cerr << "Block " << name_ << ": " << msg << std::endl << std::endl;
        exit(1);
        //return(-1);
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    const Block* find(const std::string &name) const
    //-----------------------------------------------------------------------------------
    {
        for (const auto& bl:blocks_) {
            if (bl.name_ == name) {
                return &bl;
            }
        }
        return nullptr;
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block* find(const std::string &name)
    //-----------------------------------------------------------------------------------
    {
        for (auto& bl:blocks_) {
            //std::cout << "find: " << bl.name_ << std::endl;
            if (bl.name_ == name) {
                return &bl;
            }
        }
        return nullptr;
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block& create(const std::string &name)
    //-----------------------------------------------------------------------------------
    {
        blocks_.emplace_back(name, level_+1, *this);    
        return blocks_.back();
    }
    //                                     Block

    //-----------------------------------------------------------------------------------
    Block& find_or_create(const std::string &name)
    //-----------------------------------------------------------------------------------
    {
        if (Block* block = find(name))
            return *block; 
        // Block with given name does not exist, create new
        blocks_.emplace_back(name, level_+1, *this);    
        return blocks_.back();
    }

    friend class Input;
};




//=====================================================================================
//
//                                        T A G
//
//=====================================================================================
class Tag 
{
private:
    std::vector<uint8_t> key_;
    std::vector<uint8_t> end_;
    char comm_;
    char var_;


public:
    //                                     Tag 
    //-----------------------------------------------------------------------------------
    Tag(const std::string& key="<>", const std::string& end="<end>", const char& comm='#', const char& var='$') 
        : key_(key.begin(), key.end()), end_(end.begin(), end.end()), comm_(comm), var_(var) 
    //-----------------------------------------------------------------------------------
    { 
        if (key_.empty())
            std::cout << std::endl << "!! WARNING No input keyword identifier is provided !!" << std::endl << std::endl; 
        if (key_.size() < 2)
            key_.push_back(' ');

    }

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    bool is_end(const std::vector<uint8_t>& word) const { return end_ == word ? true : false; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    bool no_end() const { return end_.empty() ? true : false; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    void set_var(const char val) { var_ = val; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    bool var_off() const { return (var_ == 0) ? true : false; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    const char comment() const { return comm_; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    const std::vector<uint8_t>& end() const { return end_; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    const uint8_t key(const int ind) const { return key_[ind]; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    const char variable() const { return var_; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    std::string all() const 
    //-----------------------------------------------------------------------------------
    { 
        return std::string(std::begin(key_), std::end(key_)) + std::string(std::begin(end_), std::end(end_)) + comm_ + var_; 
    }
    // //                                     Tag 
    // //-----------------------------------------------------------------------------------
    // friend std::ostream& operator<<(std::ostream& out, Tag& t)
    // //-----------------------------------------------------------------------------------
    // { 
    //     out << "key: " << t.key() << t.key_end() << ", end: " << t.end() << ", comment: " << t.comment() << ", var: " << t.var(); 
    //     return out; 
    // }
};

//=====================================================================================
//
//                                    I N P U T
//
//=====================================================================================
class Input
{

private:
    Tag tag_;
    Block not_found_;
    Block head_block_;
    Block * current_block_ = nullptr;
    int line_num_ = 0;
    bool recursive_ = true;
    std::string math_symbols_;

public:
    static constexpr uint8_t default_flag = 0b0000'0000; 
    static constexpr uint8_t math_off     = 0b0000'0001; 
    static constexpr uint8_t missing_ok   = 0b0000'0010;
    static constexpr uint8_t var_off      = 0b0000'0100;
    // static constexpr uint8_t skip_data    = 0b0000'1000;

    //                                     Input
    //-----------------------------------------------------------------------------------
    Input(const std::string& filename, const Tag& tag=Tag(), const uint8_t flags=default_flag)
        : tag_(tag), not_found_(), head_block_(filename, not_found_, missing_ok&flags), current_block_(&head_block_) 
    //-----------------------------------------------------------------------------------
    {         
        if (tag_.no_end())
            recursive_ = false;
        if (var_off & flags)
            tag_.set_var(0);
        if ( not (math_off & flags)) { 
            // Check for tag-conflicts with math symbols
            std::string symbols = "+-*/";
            auto tags = tag_.all(); 
            for (const auto& sym : symbols) {
                auto pos = tags.find(sym);
                if (pos == std::string::npos)
                    math_symbols_ += sym;
                else
                    std::cout << "!! WARNING Input math operation " << sym << " disabled due to conflict with tags !!" << std::endl;
            }
        }   
        parse_file();    
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    friend std::ostream& operator<<(std::ostream &out, const Input& input) { out << input.head_block_; return out; }
    //-----------------------------------------------------------------------------------

    //                                     Input
    //-----------------------------------------------------------------------------------
    const Block& operator[](const char *key) { return head_block_[key]; }
    //-----------------------------------------------------------------------------------

    //                                     Input
    //-----------------------------------------------------------------------------------
    const std::string& filename() const { return head_block_.name_; }
    //-----------------------------------------------------------------------------------



private:

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    void message(const std::string& type, const std::string& msg) 
    //-----------------------------------------------------------------------------------
    {
        std::cerr << std::endl << "*** " << type << " during reading of input-file " << filename() << " ***" << std::endl;
        std::cerr << "Line " << line_num_ << ": " << msg << std::endl << std::endl;
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    void error(const std::string& msg) 
    //-----------------------------------------------------------------------------------
    {
        message("ERROR", msg);
        exit(1);
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    template <typename T>
    void zip_flat(T a_beg, T a_end, T b_beg, T ab_beg)
    //-----------------------------------------------------------------------------------
    {    
        //std::cout << sizeof(a_beg) << std::endl;
        std::vector<int> nr(std::distance(a_beg, a_end));
        std::iota(std::begin(nr), std::end(nr), 0);    
        std::for_each(std::execution::par_unseq, std::begin(nr), std::end(nr), [&](const auto& i){
            *std::next(ab_beg, i*2)   = *std::next(a_beg, i);
            *std::next(ab_beg, i*2+1) = *std::next(b_beg, i);
        });
    }


    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    template <typename T>
    std::vector<std::vector<uint8_t>> split(const T& begin, const T& end, const std::vector<uint8_t>& delim)
    //-----------------------------------------------------------------------------------
    {
        const auto par{ std::execution::par_unseq };

        const auto size{ std::distance(begin, end) };
        std::vector<int> nr(size);
        std::iota(std::begin(nr), std::end(nr), 0);

        std::vector<int> delim_ind(size);
        std::transform(par, begin, end, std::begin(nr), std::begin(delim_ind), [&](const auto& vec, const auto& ind){
            const bool found {std::find(std::begin(delim), std::end(delim), vec) == std::end(delim)};
            return found ? -1 : ind;   
        });

        std::vector<int> word_end_ind;
        std::copy_if(std::begin(delim_ind), std::end(delim_ind), std::back_inserter(word_end_ind), [](const auto& ind){ return ind >= 0; });
        word_end_ind.push_back(size);

        std::vector<int> word_len_tmp;
        std::adjacent_difference(std::begin(word_end_ind), std::end(word_end_ind), std::back_inserter(word_len_tmp));
        std::vector<int> word_len(word_len_tmp);
        std::transform(std::begin(word_len_tmp)+1, std::end(word_len_tmp), std::begin(word_len)+1, [](const auto& len){ return len-1; });

        std::vector<std::vector<uint8_t>> words_tmp(std::size(word_len));
        std::transform(par, std::begin(word_end_ind), std::end(word_end_ind), std::begin(word_len), std::begin(words_tmp), [&begin](const auto& ind, const auto& len){
            return std::vector<uint8_t>(begin+ind-len, begin+ind);
        }); 

        std::vector<std::vector<uint8_t>> words;
        std::copy_if(std::begin(words_tmp), std::end(words_tmp), std::back_inserter(words), [](const auto& word){ 
            return std::size(word) > 0 and not std::all_of(std::begin(word), std::end(word), [](auto d){ return std::isspace(d); }); 
        });

        return words;
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    template <typename T>
    std::vector<std::vector<uint8_t>> cut(const T& begin, const T& end, uint8_t delim)
    //-----------------------------------------------------------------------------------
    {
        std::vector<std::vector<uint8_t>> tmp(std::distance(begin, end));
        std::transform(begin, end, std::begin(tmp), [&](const auto& vec)->std::vector<uint8_t>{ 
            return {std::begin(vec), std::find(std::begin(vec), std::end(vec), delim)};
        });
        // Remove empty elements
        std::vector<std::vector<uint8_t>> out;
        std::copy_if(std::begin(tmp), std::end(tmp), std::back_inserter(out), [](const auto& t){
            return std::size(t) > 0;
        });
        return out;
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    template <typename T>
    std::vector<int> matches(const T& begin, const T& end, const std::vector<uint8_t>& tokens)
    //-----------------------------------------------------------------------------------
    {
        std::vector<T> pos(std::size(tokens));
        std::transform(std::begin(tokens), std::end(tokens), std::begin(pos), [&](const auto& token){ 
            return std::find(begin, end, token);
        });
        std::vector<int> dist;
        std::transform(std::begin(pos), std::end(pos), std::back_inserter(dist), [&](const auto& p){
            return (p != end) ? std::distance(begin, p) : -1 ; 
        });
        std::vector<int> res;
        std::copy_if(std::begin(dist), std::end(dist), std::back_inserter(res), [](const auto& d){ return d>0; });
        return res;
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    template <typename T>
    int match_first(const T& begin, const T& end, const std::vector<uint8_t>& tags)
    // int match_first(const T& begin, const T& end, const U& tags)
    //-----------------------------------------------------------------------------------
    {
        std::vector<T> pos(std::size(tags));
        std::transform(std::begin(tags), std::end(tags), std::begin(pos), [&](const auto& tag){ 
            return std::find(begin, end, tag);
        });
        const auto min { *std::min_element(std::begin(pos), std::end(pos)) };
        return (min != end) ? std::distance(begin, min) : -1 ;
    }

    // //                                     Input
    // //-----------------------------------------------------------------------------------
    // //template <typename T>
    // std::vector<std::vector<uint8_t>> extract(const std::vector<uint8_t>& input, const std::vector<uint8_t>& start_tags, const std::vector<uint8_t>& stop_tags)
    // //-----------------------------------------------------------------------------------
    // {
    //     const auto start_pos { matches(std::begin(input), std::end(input), start_tags) };
    //     auto stop_pos = start_pos;
    //     std::transform(std::begin(start_pos), std::end(start_pos), std::begin(stop_pos), [&](const auto& pos){
    //         return match_first(pos, std::end(input), stop_tags);
    //     });
    //     std::vector<std::vector<uint8_t>> ext;
    //     std::transform(std::begin(start_pos), std::end(start_pos), std::begin(stop_pos), std::back_inserter(ext), [](const auto& a, const auto& b){
    //         return std::vector<uint8_t>(a,b);
    //     });
    //     return ext;
    // }

    // //                                     Input
    // //-----------------------------------------------------------------------------------
    // //template <typename T>
    // std::vector<std::vector<uint8_t>> remainder(const std::vector<uint8_t>& input, const std::vector<uint8_t>& start_tags, const std::vector<uint8_t>& stop_tags)
    // //-----------------------------------------------------------------------------------
    // {
    //     const auto start_pos { matches(std::begin(input), std::end(input), start_tags) };
    //     auto stop_pos = start_pos;
    //     std::transform(std::begin(start_pos), std::end(start_pos), std::begin(stop_pos), [&](const auto& pos){
    //         return match_first(pos, std::end(input), stop_tags);
    //     });
    //     std::vector<std::vector<uint8_t>> rem;
    //     rem.emplace_back(std::begin(input), start_pos[0]);
    //     std::transform(std::begin(start_pos)+1, std::end(start_pos), std::begin(stop_pos), std::back_inserter(rem), [](const auto& a, const auto& b){
    //         return std::vector<uint8_t>(a,b);
    //     });
    //     return rem;
    // }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------
    std::vector<uint8_t> read(const std::string& name)
    {
        std::ifstream file(name, std::ios::binary);
        if (!file) {
            std::cerr << "Error! Could not open file " + name << std::endl;
            exit(1);
        }
        std::vector<uint8_t> data((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));
        file.close();
        return data;
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    void parse_file() 
    //-----------------------------------------------------------------------------------
    {
        const auto data{ read(filename()) };    
        auto lines = split(std::begin(data), std::end(data), {'\n','\r'});
        lines = cut(std::begin(lines), std::end(lines), '#');
        for (auto& line : lines) {
            line = replace_variables(std::begin(line), std::end(line));
            const auto word{ split(std::begin(line), std::end(line), {' '}) };
            if (word[0] == tag_.end()) {
                //std::cout << "END: " << string(word[0]) << std::endl;
                current_block_ = current_block_->parent();
            } else if (word[0][0] == tag_.key(0)) {
                //std::cout << "KEY: " << string(word[0]) << std::endl;
                const auto name { string(std::begin(word[0])+1, std::end(word[0]), tag_.key(1)) };
                current_block_ = &current_block_->find_or_create(name);
            } else { 
                auto wd = std::begin(word);
                const auto name { std::string(std::begin(*wd), std::end(*wd)) };
                auto& block = current_block_->find_or_create(name);
                if (not std::atof(name.c_str()))
                    wd = std::next(wd);
                std::transform(wd, std::end(word), std::back_inserter(block.strings_), [](const auto& v){
                    return std::string(std::begin(v), std::end(v));
                });
            }
        }
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    template <typename T>
    std::vector<std::vector<uint8_t>> split(const T& begin, const T& end, const std::vector<uint8_t>& atag, const std::vector<uint8_t>& btag)
    //-----------------------------------------------------------------------------------
    {
        const auto par { std::execution::par_unseq };
        auto start_pos { matches(begin, end, atag) };
        auto stop_pos = start_pos;
        std::transform(par, std::begin(start_pos), std::end(start_pos), std::begin(stop_pos), [&](const auto& pos){
            return pos+match_first(begin+pos, end, btag);
        });
        std::vector<int> pos(1 + std::size(start_pos) + std::size(stop_pos), 0);
        zip_flat(std::begin(start_pos), std::end(start_pos), std::begin(stop_pos), std::begin(pos)+1);
        pos.push_back(std::distance(begin, end));

        std::vector<std::vector<uint8_t>> split;
        std::transform(std::begin(pos), std::end(pos)-1, std::begin(pos)+1, std::back_inserter(split), [&begin](const auto& a, const auto& b){
            return std::vector<uint8_t>(begin+a,begin+b);
        });
        return split;
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    template <typename T>
    std::vector<uint8_t> replace_variables(const T& begin, const T& end)
    //-----------------------------------------------------------------------------------
    {
        const auto line { split(begin, end, {'$'}, {'+','-','*','/',' '})};
        std::vector<std::vector<uint8_t>> repl;
        std::transform(std::begin(line), std::end(line), std::back_inserter(repl), [&](const auto& vec){
            if (vec[0]!='$')
                return vec;
            const auto* block { head_block_.find(std::string(std::begin(vec)+1,std::end(vec))) };
            const std::string var { block ? block->strings_[0] : "MISSING" };  
            return std::vector<uint8_t>(std::begin(var), std::end(var));
        });
        std::vector<uint8_t> res;
        std::for_each(std::execution::par_unseq, std::begin(repl), std::end(repl), [&](const auto& vec){ 
            res.insert(std::end(res), std::begin(vec), std::end(vec)); 
        });
        return res;
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    std::string string(const std::vector<uint8_t>& word) { return {std::begin(word), std::end(word)};}
    //-----------------------------------------------------------------------------------

    //                                     Input
    //-----------------------------------------------------------------------------------
    template <typename T>
    std::string string(const T& begin, const T& end, const uint8_t& tag) { return { begin, std::find(begin, end, tag) }; }
    //-----------------------------------------------------------------------------------

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    void read_data(std::string& word, Block& newblock) 
    //-----------------------------------------------------------------------------------
    {
        if (not tag_.var_off())
            replace_variables_with_values(word);
        if ( str_func::is_numeric(word) ) {
            newblock.push_back_number(word);
        } else {
            // String or math expression
            auto pos = word.find_first_of(math_symbols_);
            if (pos != std::string::npos) {
                // String with math symbol, check if it contains letters
                auto copy = word;
                copy.replace(pos, 1, " ");
                if (str_func::is_numeric(copy)) {
                    // Math expression
                    newblock.values_.push_back(result(copy, word[pos]));
                    return;
                }
            }
            // Plain string
            newblock.strings_.push_back(word);
        }    
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    void replace_variables_with_values(std::string& word)
    //-----------------------------------------------------------------------------------
    {
        while (1) {
            auto pos = word.find(tag_.variable());
            if (pos==std::string::npos)
                break;
            auto end = word.find_first_of("+-*/");
            if (end==std::string::npos)
                end = word.length();
            std::string var = word.substr(pos+1, end-pos-1);
            //std::cout << "var = " << var << std::endl;
            Block* found = head_block_.find(var);
            if (found) {
                word.replace(pos, end-pos, std::to_string(found->values_[0]));    
            } else {
                error("Variable " + tag_.variable() + var + " is not defined");
            }
        }
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    double result(std::string& word, char symbol)
    //-----------------------------------------------------------------------------------
    {
        std::istringstream ab(word);
        double a, b, c;
        if (ab >> a >> b) {
            switch (int(symbol)) {
                case 42: c = a*b; break;
                case 43: c = a+b; break;
                case 45: c = a-b; break;
                case 47: c = (b!=0) ? a/b : nanf(""); break;
                default: c = 0;
            }
            return c;
        } 
        std::replace(word.begin(), word.end(), ' ', symbol);
        error("The given input '" + word + "' is not supported. Math operations can not contain whitespace.");
        return 0;
    }




};


#endif /* SRC_INPUT_H_ */
