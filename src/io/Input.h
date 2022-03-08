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
#include <regex>
#include <initializer_list>


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

    //                                     Block
    // For range-based for-loops 
    //-----------------------------------------------------------------------------------
    auto begin() { return values_.begin(); }
    auto end() { return values_.end(); }
    auto begin() const { return values_.begin(); }
    auto end() const { return values_.end(); }
    auto cbegin() const { return values_.begin(); }
    auto cend() const { return values_.end(); }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    template <typename T>
    void copy_into(T& copy) const { std::copy(values_.begin(), values_.end(), std::back_inserter(copy)); }
    //-----------------------------------------------------------------------------------

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
    // template <typename T>
    // operator T() const { return T(values_.data(), values_.size()); }

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
    // Implicit conversion returning a std::vector of typename T
    // Example: std::vector<int> size = input["size-vector"];
    //
    template <typename T> 
    operator std::vector<T>() const       
    //-----------------------------------------------------------------------------------
    {
        //std::cout << "operator std::vector<T>()" << ", " << name_ << std::endl;
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
    // Only valid for vectors of size 1, i.e. a single value number
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

    // //                                     Block
    // //-----------------------------------------------------------------------------------
    // Block& operator[](const char* keyword)
    // //-----------------------------------------------------------------------------------
    // {
    //     if (blocks_.empty()) {
    //         return *this;
    //     } else {
    //         if (Block* block = find(std::string(keyword)) ) {
    //             return *block;
    //         } else {
    //             if (not missing_ok_) 
    //                 error("keyword " + std::string(keyword) + " not found!");
    //             return *not_found_;
    //         }
    //     }
    // }

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

    // //                                     Block
    // //-----------------------------------------------------------------------------------
    // std::vector<double> read_data(int size) const
    // //-----------------------------------------------------------------------------------
    // {
    //     //values_.reserve(size);
    //     std::vector<double> vec;
    //     std::ifstream file;
    //     file.open(filename_, std::ios::binary);
    //     file.seekg(pos_[0]);
    //     std::cout << name_ << ", " << pos_[0] << ", " << file.tellg() << std::endl;
    //     std::string line;
    //     std::vector<char> buf(pos_[1]-pos_[0]);
    //     file.read(&buf[0], buf.size());
    //     std::string data(buf.begin(), buf.end());
    //     std::cout << data << std::endl;
    //     std::istringstream stream(data);
    //     double value;
    //     while (stream >> value)
    //         vec.push_back(value);
    //     return vec;
    // }

    // //                                     Block
    // //-----------------------------------------------------------------------------------
    // const std::vector<double>& read_binary_data()
    // //-----------------------------------------------------------------------------------
    // {
    //     std::ifstream file;
    //     file.open(parent_->name_, std::ios::binary);
    //     int size = pos_[1] - pos_[0];
    //     values_.reserve(size);
    //     file.seekg(pos_[0]);
    //     //file.read(&data[0], size);
    //     return values_;
    // }

    //                                     Block
    //-----------------------------------------------------------------------------------
    friend std::ostream& operator<<(std::ostream& out, const Block& block)
    //-----------------------------------------------------------------------------------
    {
        std::string indent = std::string(block.level_*3, ' ');
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
    std::string key_;
    std::string end_;
    char comm_ = '#';
    char var_ = '$';


public:
    
    //                                     Tag 
    //-----------------------------------------------------------------------------------
    Tag() : key_("<>"), end_("<end>") { }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    Tag(const std::initializer_list<std::pair<std::string, std::string>>& ilist)
    //-----------------------------------------------------------------------------------
    {
        for (const auto& [name, val] : ilist) {
            if (name == "keyword")
                key_ = val;
            else if (name == "end")
                end_ = val;
            else if (name == "comment")
                comm_ = val[0];
            else if (name == "variable")
                var_ = val[0];
            else {
                std::cerr << "ERROR Unknown Input-tag: " << val << std::endl;
                exit(1);
            }
        }
    }

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    char key() const { return (key_.empty()) ? 0 : key_[0]; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    char key_end() const { return (key_.length()==2) ? key_[1] : 0; }    
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    const std::string& end() const { return end_; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    bool is_end(const std::string& word) const 
    //-----------------------------------------------------------------------------------
    {
        if (end_.empty())
            return false;
        if (end_.compare(word) == 0) 
            return true;
        else
            return false;
    }

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    // const std::string& comment() const { return map_.at("comment"); }
    const char comment() const { return comm_; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    // const std::string& var() const { return map_.at("variable"); }
    const char var() const { return var_; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    void set_var(const char val) { var_ = val; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    bool var_off() const { return (var_ == 0) ? true : false; }
    //-----------------------------------------------------------------------------------

    // //                                     Tag 
    // //-----------------------------------------------------------------------------------
    // const std::string& operator[](const std::string& var) { return map_[var]; }
    // //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    std::string all() const { return key_ + end_ + var_ + comm_; }
    //-----------------------------------------------------------------------------------
    
    //                                     Tag 
    //-----------------------------------------------------------------------------------
    friend std::ostream& operator<<(std::ostream& out, Tag& t)
    //-----------------------------------------------------------------------------------
    { 
        out << "key: " << t.key() << t.key_end() << ", end: " << t.end() << ", comment: " << t.comment() << ", var: " << t.var(); 
        return out; 
    }
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
    std::ifstream infile_;

public:
    static constexpr uint8_t default_flag = 0b0000'0000; 
    static constexpr uint8_t math_off     = 0b0000'0001; 
    static constexpr uint8_t missing_ok   = 0b0000'0010;
    static constexpr uint8_t var_off      = 0b0000'0100;

    //                                     Input
    //-----------------------------------------------------------------------------------
    Input(const std::string& filename, const Tag& tag=Tag(), const uint8_t flags=default_flag)
        : tag_(tag), not_found_(), head_block_(filename, not_found_, missing_ok&flags), current_block_(&head_block_) 
    //-----------------------------------------------------------------------------------
    {         
        if (tag_.end().empty())
            recursive_ = false;
        if (var_off & flags)
            tag_.var_off();
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
        process_input();    
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
    int file_pos() { return infile_.tellg(); }
    //-----------------------------------------------------------------------------------

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    void process_input() 
    //-----------------------------------------------------------------------------------
    {
        open();
        std::string word;
        std::string line;
        while ( std::getline(infile_, line) ) {
            //remove_space(line);
            remove_comments(line);
            if (line.empty()) {
                line_num_++;
                continue;
            }
            // std::cout << "line: |" << line << "|" << std::endl;            
            std::istringstream iss_line(line);
            if ( iss_line >> word ) { 
                line_num_++;              
                if (tag_.is_end(word)) {
                    // std::cout << "END: " << word << std::endl;
                    current_block_ = current_block_->parent_;
                } else if (is_keyword(word)) {
                    if (not recursive_ and current_block_->parent_) {
                        current_block_ = current_block_->parent_;
                    }
                    // std::cout << "KEY: " << word << std::endl;
                    read_keyword(word, iss_line);
                } else {
                    // we are inside a block or sub-block
                    // std::cout << "CONTENT: " << word << ", " << iss_line.str() << std::endl;
                    read_block_content(word, iss_line);
                }
            }
        }
        infile_.close();
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    bool is_keyword(const std::string &word) 
    //-----------------------------------------------------------------------------------
    {
        size_t a = word.find(tag_.key());
        if (a != std::string::npos) {
            // if (key_tag_.length()>1) {
            //     size_t b = word.find(end);
            // if (b < std::string::npos) {
            //     if (a>0 || b<word.length()-1) {
            //         error("Keyword identifiers inside keyword " + word);
            //     }
            // }
            return true;
        }
        return false;
    }



    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    void read_keyword(std::string& word, std::istringstream& stream) 
    //-----------------------------------------------------------------------------------
    {
        auto name = remove_key_tags(word);
        if (name.empty())
            return;
        Block& newblock = current_block_->find_or_create(name);  // find the matching state
        //current_block_ = &newblock; 
        bool outside = false;
        Block* data_block = nullptr;
        bool end_passed = false;
        while (stream >> word) {
            if (not tag_.key_end()) {
                // No end-tag given, add word to the block
                read_data(word, newblock);
            } else {
                if (tag_.is_end(word)) {
                    return;
                }
                if (end_passed) {
                    data_block = &newblock.create(word);
                    outside = true;
                    end_passed = false;
                }
                if (outside) { 
                    read_data(word, *data_block);
                } else {
                    newblock.datatype_ = word;                
                    // Check for end-tag 
                    auto pos = word.find(tag_.key_end());
                    if (pos != std::string::npos) {
                        // We at the end of the keyword 
                        end_passed = true;
                    }
                }    
            }
        }
        current_block_ = &newblock; 
    }


    //                                     Input
    //-----------------------------------------------------------------------------------
    // If the first word is a number we assume it is a data-value
    // and use the line-number of the block to name the block
    //
    //-----------------------------------------------------------------------------------
    void read_block_content(std::string& word, std::istringstream& stream) 
    {
        // if (skip_data_)
        //     return;
        Block& newblock = current_block_->find_or_create(word);

        // Remove first word from stream if it is a string (used to name the block)
        if (not str_func::is_numeric(word)) 
            stream >> word;

        // Process the rest of the line
        while ( 1 ) {
            read_data(word, newblock);
            if ( not (stream >> word) )
                break;
        }
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    void read_data(std::string& word, Block& newblock) 
    //-----------------------------------------------------------------------------------
    {
        if (not tag_.var_off())
            replace_variables_with_values(word);
        if ( str_func::is_numeric(word) ) {
            // Numeric value, add it
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
            auto pos = word.find(tag_.var());
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
                error("Variable " + std::string(1, tag_.var()) + var + " is not defined");
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



    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------
    void open() //const std::string& filename, std::ifstream& file)
    {
        infile_.open(filename(), std::ios::binary);
        if (not infile_) {
            std::cerr << "Error! Could not open file " + filename() << std::endl;
            exit(1);
        }
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------
    void remove_space(std::string &str) 
    {
        // remove leading and trailing space
        size_t a = str.find_first_not_of(' ');
        if (a == std::string::npos) a = 0;
        str = str.substr(a, str.find_last_not_of(' ') + 1 - a);
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------
    void remove_comments(std::string &str) 
    {
        // only keep line up to comment character
        size_t a = str.find_first_of(tag_.comment());
        if (a < std::string::npos)
            str = str.substr(0,a);

        // // remove /* */ comments
        // a = str.find("/*");
        // if (a < std::string::npos) {
        //     size_t b = str.find("*/");
        //     if (b < std::string::npos)
        //         str.erase(a, b-a+2);
        //     else
        //         std::cerr << "Warning! Missing end-comment (*/) in input-file" << std::endl;
        // }
    }


    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    const std::string& remove_key_tags(std::string& name)
    //-----------------------------------------------------------------------------------
    {
        size_t a = name.find(tag_.key());
        if (a != std::string::npos) {
            name.erase(name.begin()+a);
        }
        // Remove end tag if given
        // if (key_tag_.length() > 1) {
        if (tag_.key_end()) {
            size_t b = name.find(tag_.key_end());
            if (b != std::string::npos) {
                name.erase(name.begin()+b, name.end());
            }
        }
        return name;
    }

};


#endif /* SRC_INPUT_H_ */
