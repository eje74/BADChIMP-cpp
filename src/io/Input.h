/*
 * Input.h
 *
 *  Created on: 4. mar. 2015
 *      Author: janlv
 */

#ifndef SRC_INPUT_H_
#define SRC_INPUT_H_

#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <stack>
#include <queue>
//#include <deque>
#include <typeinfo>
#include <numeric>
#include <algorithm>
#include <memory>
#include "../lbsolver/LBlatticetypes.h"


//------------------------------------------------------------
// Demonstration of the Input-class
//
// Contents of file example.file
//
// # use # to comment whole lines or part of a line
// set version 1.1 # one-liner set commands
// <math>
//    pi 3.14
//    e  2.72
//    <series>
//       fibo 1 1 2 3 5 8
//       <prime int>   # unless 'int' or 'char' is specified, a 'double' datatype is assumed
//          2 3 5 7
//          11 13 17 19
//       <end>
//    <end>
// <end>
// <binary char>    # use a char datatype to read digit by digit
//    0101
//    1110
// <end>
//
// How to access contents of example.file in code:
//
// Input input("example.file");
// double version = input["version"];                      // version = 1.1
// double pi = input["math"]["pi"];                        // pi = 3.12
// std::vector<int> fibo = input["math"]["series"]["fibo"]; // fibo = 2, 3, 5, 7, 11, 13, 17, 19
// std::vector<int> binary = input["binary"];              // binary = 0, 1, 0, 1, 1, 1, 1, 0
//


namespace my_string
{
    //-----------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------
    bool is_numeric(const std::string& str)
    {
        std::istringstream ss(str);
        double dbl;
        ss >> dbl;      // try to read the number

        //std::cout << "in: (" << str << ")"; 
        if (!ss.fail() && ss.eof()) {
            //std::cout << " is NUMBER" << std::endl;
            return true;  // is-a-number
        } else {
            //std::cout << " is NOT A NUMBER" << std::endl;
            return false; // not-a-number
        }
    }

    //-----------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------
    bool is_string(const std::string& str)
    {
        return !is_numeric(str);
    }

    // //-----------------------------------------------------------------------------------
    // //
    // //-----------------------------------------------------------------------------------
    // const std::string increment(std::string &s) { return std::to_string(std::stoi(s)+1); }
}

//=====================================================================================
//
//                                    B L O C K 
//
// A block is the content between an opening and a closing tag
// A block can contain blocks recursively
// A block is accessed by the []-operator
//=====================================================================================

//template <typename T>
class Block
{
public:
    int level_ = 0;                 // used for indenting output
    std::string name_;              // name of the block given by the first string on the line. If no string, the line number is used
    std::vector<Block> blocks_;     // a block is recursive and can contain other (sub-) blocks
    std::vector<double> values_;    // the values of the block, could be a template
    std::vector<std::string> strings_;    // the values of the block, could be a template
    std::string datatype_;          // name of the datatype, probably obsolete if a template class is used
    Block* parent_ = nullptr;

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Constructors
    //
    Block() : level_(0), name_(), blocks_(), values_(), strings_(), datatype_(), parent_(nullptr) {}
    Block(const std::string &name, int level, Block* parent) 
        : level_(level), name_(name), blocks_(), values_(), strings_(), datatype_(), parent_(parent)
    //-----------------------------------------------------------------------------------
    {
        if ( my_string::is_numeric(name_) ) {
            // use line-number as name
            name_ = "0";
            if (parent != nullptr)
                name_ = std::to_string(parent->blocks_.size());
        }
    }


    //                                     Block
    //-----------------------------------------------------------------------------------
    void info() const
    //-----------------------------------------------------------------------------------
    {
        std::cout << name_ << ", level = " << level_ << ", blocks = ";
        for (const auto& b : blocks_)
            std::cout << b.name_ << ", ";
        std::cout << "values = " << values_.size();
        std::cout << ", strings = " << strings_.size() << ", datatype = " << datatype_ << ", parent = " << ((parent_==nullptr)? ("nullptr") : parent_->name_);
        std::cout << std::endl;                
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    operator int() const { return static_cast<int>(values_[0]); }
    //operator int() const { return std::get<int>(values_[0]); }
    operator double() const { return values_[0]; }
    //operator T() const { return values_[0]; }
    operator std::string() const { return strings_[0]; }
    // operator int() const { std::cout << name_ << " Input::int()" << std::endl; return static_cast<int>(values_[0]); }
    //-----------------------------------------------------------------------------------
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    // Block<T>& at() 
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
    int nrows_match(const std::string &pattern)
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
    int nrows_not_match(const std::string &pattern) { return nrows()-nrows_match(pattern); }
    //-----------------------------------------------------------------------------------
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    // Implicit conversion returning a std::vector of typename T
    // Example: std::vector<int> size = input["size-vector"];
    //
    template <typename T> 
    operator std::vector<T>()            
    //-----------------------------------------------------------------------------------
    {
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
            // only one line, return vector
            return(std::vector<T>(values_.begin(), values_.end()));
        }
    }   

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Implicit conversion returning a std::vector of std::string
    //
    operator std::vector<std::string>() { return(std::vector<std::string>(strings_.begin(), strings_.end())); }   
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Implicit conversion returning a std::valarray of typename T
    // Example: std::valarray<int> size = input["size-vector"];
    //
    template <typename T> 
    operator std::slice_array<T>()            
    //-----------------------------------------------------------------------------------
    {
        std::vector<T> tmp = *this;
        std::valarray<T> varr(tmp.data(), tmp.size());
        return(varr[std::slice(0, varr.size(), 1)]);
    }

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
    //std::vector<double> column(const int n)    
    template <typename T>
    std::vector<T> column(const int n)    
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
    std::vector< std::vector<T> > matrix()            
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
    std::vector< std::vector<T> > symmetric_matrix()
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
    std::vector<std::string> names(void)       
    //-----------------------------------------------------------------------------------
    {
        //std::vector<std::string> *var_names = new std::vector<std::string>();
        std::vector<std::string> var_names;// = new std::vector<std::string>();
        for (size_t i=0; i<blocks_.size(); ++i) {
            var_names.emplace_back(blocks_[i].name_);
        }
        return var_names;
    }
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    Block & operator[](const char *keyword)
    //-----------------------------------------------------------------------------------
    {
        if (blocks_.empty()) {
            return *this;
        } else {
            Block* ret = find(std::string(keyword));
            if (ret==nullptr) {
                std::cerr << "ERROR in Block::operator[], keyword " << keyword << " not found!" << std::endl;
                exit(1);
            }
            return *ret;
        }
    }


    //                                     Block
    //-----------------------------------------------------------------------------------
    //double operator[](const int ind) {
    //template <typename T>
    double operator[](const int ind) {
    //-----------------------------------------------------------------------------------
        if (int(values_.size())<ind+1) {
            std::cerr << std::endl << "ERROR reading input-file: '" << name_ << "': index outside limit, " << ind << " > " << values_.size()-1 << std::endl << std::endl;
            std::cerr << values_[0] << std::endl;
            exit(1);
        }
        return values_[ind];
    }


    //                                     Block
    //-----------------------------------------------------------------------------------
    friend std::ostream& operator<<(std::ostream& out, const Block& block)
    //-----------------------------------------------------------------------------------
    {
        std::string indent = std::string(block.level_*3, ' ');
        out << indent+block.name_ << ":";
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
    Block* find(const std::string &name)
    //-----------------------------------------------------------------------------------
    {
        for (auto& bl:blocks_) {
            if (bl.name_ == name) {
                return &bl;
            }
        }
        return nullptr;
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block& find_or_create(const std::string &name)
    //-----------------------------------------------------------------------------------
    {
        Block* block = find(name);
        if (block != nullptr) 
            return *block;
        // Block with given name does not exist, create new
        blocks_.emplace_back(name, level_+1, this);    
        return blocks_.back();
    }

    // //                                     Block
    // //-----------------------------------------------------------------------------------
    // template <typename T>
    // void print_vec(std::vector<T> &vec)
    // //-----------------------------------------------------------------------------------
    // {
    //     for (const auto& v:vec)
    //         std::cout << v << " ";
    //     std::cout << std::endl;
    // }

    friend class Input;
};


//template <typename T>
//=====================================================================================
//
//                                    I N P U T
//
//=====================================================================================
class Input
{
private:
    std::string key_start_id_ = "<";
    std::string key_end_id_ = ">";
    std::string end_word_ = "end";
    std::string set_word_ = "set";
    std::ifstream infile_;
    Block head_block_;
    std::queue<std::istringstream*> input_;
    Block * current_block_ = nullptr;

public:
    //                                     Input
    //-----------------------------------------------------------------------------------
    Input(const std::string filename) : infile_(), head_block_(), input_() 
    //-----------------------------------------------------------------------------------
    {
        end_word_ = key_start_id_ + end_word_ + key_end_id_;
        current_block_ = &head_block_;
        set_input_from_file(filename);
        process_input();
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    friend std::ostream& operator<<(std::ostream &out, const Input& input) 
    //-----------------------------------------------------------------------------------
    { 
        out << input.head_block_;
        return out;
    }

    const Block& head_block() const {return head_block_;}

    //                                     Input
    //-----------------------------------------------------------------------------------
    Block& operator[](const char *key)
    //-----------------------------------------------------------------------------------
    {
        std::string keyword(key);
        Block* ret = head_block_.find(keyword);
        if (ret == nullptr) {
            std::cout << "Error! keyword " << keyword << " not found!" << std::endl;
            exit(1);
        } else {
            return *ret;
        }
    }


private:
    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    void process_input() 
    //-----------------------------------------------------------------------------------
    {
        std::string word;
        std::istringstream *line;
        while ( !input_.empty() ) {
            line = input_.front();
            //std::cout << "line: (" << line->str() << ")" << std::endl;            
            if ( (*line) >> word ) {                
                if (word == end_word_) {
                    // std::cout << "END: " << word << std::endl;
                    current_block_ = current_block_->parent_;
                } else if (is_keyword(word)) {
                    // std::cout << "KEY: " << word << std::endl;
                    read_keyword(word, line);
                } else {
                    // we are inside a block or sub-block
                    // std::cout << "CONTENT: " << word << std::endl;
                    read_block_content(word, line);
                }
            }
            input_.pop();
        }
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    // Fill queue of stringstreams from input-file
    //
    void set_input_from_file(const std::string &filename)
    //-----------------------------------------------------------------------------------
    {
        //std::cout << "Input::init: using file " << filename << std::endl;
        head_block_.name_ = filename;
        open(filename);
        std::string line;
        while ( std::getline(infile_, line) ) {
            //std::cout << "A: (" << line << ")" << std::endl;
            remove_space(line);
            remove_comments(line);
            //std::cout << "B: (" << line << ")" << std::endl;
            if (line.empty())
            continue;
            //std::cout << '(' << line.c_str() << ')' << std::endl;
            //std::cout << '(' << ')'; // << std::endl;
            //printf("(%s)\n",line.c_str());

            input_.push(new std::istringstream(line));
        }
        infile_.close();
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    void read_set(std::istringstream *stream) 
    //-----------------------------------------------------------------------------------
    {
        std::string val;
        std::string var;
        // Allow number*number syntax
        // std::string var, val2;
        // (*stream) >> var >> val >> val2;
        // //std::cout << "var: " << var << ", val: " << val << ", val2: " << val2 << std::endl;
        // if (val2.size()>0) {
        //   if (val2[0]=='*') {
        //     val2.erase(0,1);
        //     val *= std::stod(val2);
        //   }
        //   //std::cout << "VAL: " << val << std::endl;
        // }
        (*stream) >> var >> val;
        Block& block = current_block_->find_or_create(var);
        if (my_string::is_numeric(val)) {
            block.values_.push_back(std::stod(val));
        } else {
            block.strings_.push_back(val);
        }
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    void read_keyword(std::string &word, std::istringstream *stream) 
    //-----------------------------------------------------------------------------------
    {
        Block& newblock = current_block_->find_or_create(remove_key_tags(word));  // find the matching state
        while ((*stream) >> word) {
            newblock.datatype_ = remove_key_tags(word);
        }
        current_block_ = &newblock; 
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    // If the first word is a number we assume it is a data-value
    // and use the line-number of the block to name the block
    //
    //-----------------------------------------------------------------------------------
    void read_block_content(std::string &word, std::istringstream *stream) 
    {
        Block *parent = current_block_;
        Block& newblock = parent->find_or_create(word);

        if (my_string::is_string(word)) {
            if ( !((*stream) >> word) ) {
                std::cerr << "ERROR in Input: Missing values for " << word << std::endl;
                exit(1);
            }   
        }

        // process the rest of the line
        while ( 1 ) {
            if ( my_string::is_numeric(word) ) {
                if (parent->datatype_ == "char") {
                    push_back_word<char>(word, newblock);
                } else if (parent->datatype_ == "int") {
                    push_back_word<int>(word, newblock);
                } else {
                    push_back_word<double>(word, newblock);
                }
            } else {
                newblock.strings_.push_back(word);
                // create new blocks for each additional word on a line
                //newblock = newblock.find_or_create(word);
            }
            if ( !((*stream) >> word) )
                break;
        }
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------
    bool is_keyword(const std::string &word) 
    {
        size_t a = word.find_first_of(key_start_id_);
        if (a < std::string::npos) {
            size_t b = word.find_first_of(key_end_id_);
            //std::cout << "a,b:" << a << "," << b << std::endl;
            if (b < std::string::npos) {
                if (a>0 || b<word.length()-1) {
                    std::cerr << "Input error! Keyword identifiers found inside keyword: "+word << std::endl;
                    exit(-1);
                }
            //return true;
            }
            return true;
            //    else {
            //      std::cerr << "Input error! Missing " + key_end_id + " in keyword: "+word << std::endl;
            //      exit(-1);
            //    }
        }
        return false;
    }



    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------
    int open(const std::string& filename)
    {
        infile_.open(filename.c_str());
        if (!infile_) {
            std::cerr << "Error! Could not open file " + filename << std::endl;
            exit(1);
        }
        return 1;
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
        // only keep line up to '#' character
        size_t a = str.find_first_of('#');
        if (a < std::string::npos)
            str = str.substr(0,a);

        // remove /* */ comments
        a = str.find("/*");
        if (a < std::string::npos) {
            size_t b = str.find("*/");
            if (b < std::string::npos)
                str.erase(a, b-a+2);
            else
                std::cerr << "Warning! Missing end-comment (*/) in input-file" << std::endl;
        }
    }


    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------
    const std::string& remove_key_tags(std::string& name)
    {
        size_t a = name.find_first_of(key_start_id_);
        if (a < std::string::npos) {
            name.erase(name.begin()+a);
        }
        size_t b = name.find_first_of(key_end_id_);
        if (b < std::string::npos) {
            name.erase(name.begin()+b, name.end());
        }
        return name;
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------
    template <typename T>
    // void push_back_word(const std::string& word, Block* newblock) 
    void push_back_word(const std::string& word, Block& newblock) 
    {
        T w;
        std::istringstream iss(word);
        while (iss >> w) {
            if (typeid(w) == typeid(char))
                w -= '0';
            // newblock->values_.push_back(w);
            newblock.values_.push_back(w);
        }
    }
};


#endif /* SRC_INPUT_H_ */
