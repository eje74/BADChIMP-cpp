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
#include <typeinfo>
#include <numeric>
#include <algorithm>
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
// std::vector<int> fibo = input["math"]["series"]["fibo"] // fibo = 2, 3, 5, 7, 11, 13, 17, 19
// std::vector<int> binary = input["binary"];              // binary = 0, 1, 0, 1, 1, 1, 1, 0
//


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
    int level_ = 0;                 // used for indenting output
    std::string name_;              // name of the block given by the first string on the line. If no string, the line number is used
    std::vector<Block*> blocks_;    // a block is recursive and can contain other (sub-) blocks
    std::vector<double> values_;    // the values of the block, could be a template
    std::string datatype_;          // name of the datatype, probably obsolete if a template class is used


    //                                     Block
    //-----------------------------------------------------------------------------------
    // Constructors
    //
    Block() {}
    Block(const std::string &name, int level) : level_(level), name_(name) {}
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Destructor
    //
    ~Block() { for (auto b:blocks_) delete b; }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    operator int() const { return static_cast<int>(values_[0]); }
    operator double() const { return values_[0]; }
    operator std::string() const { return std::to_string(values_[0]); }
    // operator int() const { std::cout << name_ << " Input::int()" << std::endl; return static_cast<int>(values_[0]); }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Return number of rows
    size_t nrows(void) const { return blocks_.size(); }
    //-----------------------------------------------------------------------------------
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    // Return number of columns
    size_t ncols(void) const { return values_.size(); }
    //-----------------------------------------------------------------------------------
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    // return number of rows matching given string
    int nrows_match(const std::string &pattern)
    //-----------------------------------------------------------------------------------
    {
        int m=0;
        for (const auto& b:blocks_) {
            if (b->name_ == pattern)
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
        if (blocks_.size()>0) {
            // If the block contains several unnamed lines, return a flattened
            // vector of all lines of the block.
            // See <prime> or <binary> in the example.file on top.
            std::vector<T> vec;
            vec.reserve(nrows()*blocks_[0]->ncols());
            for (const auto& bl:blocks_) {
                vec.insert(vec.end(), bl->values_.begin(), bl->values_.end());
            }
            return vec;
        } else {
            // only one line, return vector
            return(std::vector<T>(values_.begin(), values_.end()));
        }
    }   

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
    std::vector<double> get_row() const { return(values_); }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Create and return a 1D VECTOR of column n of a block
    //
    std::vector<double> get_column(const int n)    
    //-----------------------------------------------------------------------------------
    {
        std::vector<double> col;
        for (const auto& bl:blocks_) {
            if ( (n+1)>int(bl->values_.size())) {
            std::cerr << "ERROR in Block::get_column: index " << n
                << " beyond limit " << bl->values_.size()-1 << std::endl;
            exit(1);
            }
            col.push_back(bl->values_[n]);
        }
        return col;
    }

    
    //                                     Block
    //-----------------------------------------------------------------------------------
    // Create and return a MATRIX (2D vector) of size nrows by ncols
    //
    std::vector< std::vector<double> > get_2D_matrix()            
    //-----------------------------------------------------------------------------------
    {
        std::vector< std::vector<double> > mat(nrows(), std::vector<double>(blocks_[0]->ncols()));
        int n = 0;
        for (auto &row : mat) {
            row = blocks_[n++]->get_row();
        }
        return mat;
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Create and return a SYMMETRIC MATRIX of size nrows by ncols
    // where the lower triangle is overwritten by the upper triangle
    // 
    std::vector< std::vector <double> > get_symmetric_2D_matrix()
    //-----------------------------------------------------------------------------------
    {
        std::vector< std::vector<double> > mat = get_2D_matrix();
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
    std::vector<std::string> & get_names(void)       
    //-----------------------------------------------------------------------------------
    {
        std::vector<std::string> *var_names = new std::vector<std::string>();
        for (int i=0; i<(int)blocks_.size(); i++) {
            var_names->push_back(blocks_[i]->name_);
        }
        return (*var_names);
    }
    
    //                                     Block
    //-----------------------------------------------------------------------------------
    Block & operator[](const std::string &keyword)
    //-----------------------------------------------------------------------------------
    {
        if (blocks_.empty()) {
            return *this;
        } else {
            Block *ret = find(keyword);
            if (ret==nullptr) {
                std::cerr << "ERROR in Block::operator[], keyword " << keyword << " not found!" << std::endl;
                exit(1);
            }
            return *ret;
        }
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block & operator[](const char *key) { return (*this)[std::string(key)]; }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
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
    //std::ostream& operator<<(std::ostream& out){ return out << values[0];};
    void print(void) 
    //-----------------------------------------------------------------------------------
    {
        std::string indent = std::string(level_*3, ' ');
        std::cout << indent+name_ << ":";
        for (const auto& val : values_) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // recursive call
        for (const auto& block : blocks_) {
            block->print();
        }
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    // Print value array as a space separated list
    //-----------------------------------------------------------------------------------
    void print_value(const std::string &keyword)
    {
        for (const auto& val : values_) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }


private:
    // void remove_key_identifyer(std::string &start_id, std::string &end_id);

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block * find(const std::string &name)
    //-----------------------------------------------------------------------------------
    {
        //std::cout << "Trying to find " << key << " in block " << name << ".....";
        //std::vector<Block<T>*>::iterator it;
        //for(it=blocks.begin(); it!=blocks.end(); it++) {
        for (const auto& b:blocks_) {
            if (b->name_ == name) {
                //std::cout << "success!!" << std::endl;
                return b;
            }
        }
        return nullptr;
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block * find_or_create(const std::string &name)
    //-----------------------------------------------------------------------------------
    {
        Block *b = find(name);
        if (b != NULL) return b;
        // did not exist, create new
        //std::cout << "Creating new block with name " << name << std::endl;
        blocks_.push_back(new Block(name, level_+1));
        return blocks_.back();
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    template <typename T>
    void print_vec(std::vector<T> &vec)
    //-----------------------------------------------------------------------------------
    {
        for (const auto& v:vec)
            std::cout << v << " ";
        std::cout << std::endl;
    }

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
    std::string key_start_id = "<";
    std::string key_end_id = ">";
    std::string end_word = "end";
    std::string set_word = "set";
    std::ifstream infile;
    std::stack<Block*> current_block;
    Block *head_block;
    std::queue<std::istringstream*> input;

public:
    //                                     Input
    //-----------------------------------------------------------------------------------
    Input(const std::string filename) : infile(), current_block(), head_block(new Block()), input() 
    //-----------------------------------------------------------------------------------
    {
        end_word = key_start_id + end_word + key_end_id;
        current_block.push(head_block);
        set_input_from_file(filename);
        process_input();
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    ~Input() 
    //-----------------------------------------------------------------------------------
    {
        current_block.pop();
        delete head_block;
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    void print(void)
    //-----------------------------------------------------------------------------------
    {
        head_block->print();
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    Block & operator[](const char *key)
    //-----------------------------------------------------------------------------------
    {
        Block *ret;
        std::string keyword(key);
        //keyword.assign(key);
        //std::cout << "Searching for keyword " << keyword << std::endl;
        ret = head_block->find(keyword);
        if (ret == NULL) {
            std::cout << "Error! keyword " << keyword << " not found!" << std::endl;
            exit(1);
        } else {
            return(*ret);
        }
    }


private:
    //                                     Input
    //-----------------------------------------------------------------------------------
    void process_input() 
    //-----------------------------------------------------------------------------------
    {
        std::string word;
        std::istringstream *line;
        while ( !input.empty() ) {
            line = input.front();
            //std::cout << "line: (" << line->str() << ")" << std::endl;            
            if ( (*line) >> word ) {                
                if (word == set_word) {
                    //std::cout << "SET: " << word << std::endl;
                    read_set(line);
                } else if (word == end_word) {
                    //std::cout << "END: " << word << std::endl;
                    current_block.pop();
                } else if (is_keyword(word)) {
                    //std::cout << "KEY: " << word << std::endl;
                    read_keyword(word, line);
                } else {
                    // we are inside a block or sub-block
                    //std::cout << "CONTENT: " << word << std::endl;
                    read_block_content(word, line);
                }
            }
            input.pop();
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
        head_block->name_ = filename;
        open(filename);
        std::string line;
        while ( std::getline(infile, line) ) {
            //std::cout << "A: (" << line << ")" << std::endl;
            remove_space(line);
            remove_comments(line);
            //std::cout << "B: (" << line << ")" << std::endl;
            if (line.empty())
            continue;
            //std::cout << '(' << line.c_str() << ')' << std::endl;
            //std::cout << '(' << ')'; // << std::endl;
            //printf("(%s)\n",line.c_str());

            input.push(new std::istringstream(line));
        }
        infile.close();
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    void read_set(std::istringstream *stream) 
    //-----------------------------------------------------------------------------------
    {
        double val;
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
        Block *block = current_block.top()->find_or_create(var);
        block->values_.push_back(val);
    }

    void read_keyword(std::string &word, std::istringstream *stream) 
    {
        //size_t a = word.find_first_of(key_start_id);
        //std::cout << "a: " << a << std::endl;
        Block *newblock = current_block.top()->find_or_create(remove_key_tags(word));  // find the matching state
        current_block.push(newblock);        // add to stack
        while ((*stream) >> word) {
            newblock->datatype_ = remove_key_tags(word);
            //std::cout << newblock->datatype << std::endl;
            //newblock->values.push_back(std::stof(word));
            //    if ( is_numeric(word) ) {    // if number: add to range (int)
            //      newblock->range.push_back(std::stoi(word)); // std::stoi is C++11
            //    } else {
            //      std::cout << "WARNING! Input ignores non-numeric words after keywords: "
            //          << "'" << word << "'" << " given after " << current_block.top()->name <<std::endl;
            //    }
        }
    }
    void read_end(void);

    //---------------------
    // If the first word is a number we assume it is a data-value
    // and use the line-number of the block to name the block
    //
    //---------------------
    void read_block_content(std::string &word, std::istringstream *stream) 
    {
        Block *parent = current_block.top();
        
        //if (parent->word_length)
        //  std::cout << parent->name << ", word_length" << parent->word_length << std::endl;
        // create name for the new block
        std::string name;
        if ( is_numeric(word) ) {
            // use line-number as name
            name = "0";
            if (parent->blocks_.size()>0)
                name = inc_string(parent->blocks_.back()->name_);
        } else {
            // use given word as name
            name = word;
            //newblock = current_block.top()->find_or_create(word);
            if ( !((*stream) >> word) ) {
                std::cerr << "ERROR in Input: Missing values for " << name << std::endl;
                exit(1);
            } 
        }
        
        // create new block
        Block *newblock = parent->find_or_create(name);

        // process the rest of the line
        while ( 1 ) {
            if ( is_numeric(word) ) {
                if (parent->datatype_ == "char") {
                    push_back_word<char>(word, newblock);
                } else if (parent->datatype_ == "int") {
                    push_back_word<int>(word, newblock);
                } else {
                    push_back_word<double>(word, newblock);
                }
            } else {
                // create new blocks for each additional word on a line
                newblock = newblock->find_or_create(word);
            }
            if ( !((*stream) >> word) )
                break;
        }
    }

    bool is_keyword(const std::string &word) 
    {
        size_t a = word.find_first_of(key_start_id);
        if (a < std::string::npos) {
            size_t b = word.find_first_of(key_end_id);
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


    bool is_numeric(const std::string& str) const
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

    int open(const std::string& filename)
    {
        infile.open(filename.c_str());
        if (!infile) {
            std::cerr << "Error! Could not open file " + filename << std::endl;
            exit(1);
        }
        return 1;
    }

    void remove_space(std::string &str) 
    {
        // remove leading and trailing space
        size_t a = str.find_first_not_of(' ');
        if (a == std::string::npos) a = 0;
        str = str.substr(a, str.find_last_not_of(' ') + 1 - a);
    }

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

    const std::string inc_string(std::string &s) const 
    {
        return std::to_string(std::stoi(s)+1);
    }

    //                                     Input
    //-----------------------------------------------------------------------------------
    const std::string& remove_key_tags(std::string& name)
    //-----------------------------------------------------------------------------------
    {
        size_t a = name.find_first_of(key_start_id);
        if (a < std::string::npos) {
            name.erase(name.begin()+a);
        }
        size_t b = name.find_first_of(key_end_id);
        if (b < std::string::npos) {
            name.erase(name.begin()+b, name.end());
        }
        return name;
    }

    template <typename T>
    void push_back_word(const std::string& word, Block* newblock) 
    {
        T w;
        std::istringstream iss(word);
        while (iss >> w) {
            if (typeid(w) == typeid(char))
                w -= '0';
            newblock->values_.push_back(w);
        }
    }
};


#endif /* SRC_INPUT_H_ */
