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
//------------------------------------------------------------


//------------------------------------------------------------
// A block is the content between an opening and a closing tag
// A block can contain blocks recursively
// A block is accessed by the []-operator
//------------------------------------------------------------
class Block {
 public:
  int level_ = 0;                 // used for indenting output
  std::string name_;              // name of the block given by the first string on the line. If no string, the line number is used
  std::vector<Block*> blocks_;    // a block is recursive and can contain other (sub-) blocks
  std::vector<double> values_;    // the values of the block, could be a template
  std::string datatype_;          // name of the datatype, probably obsolete if a template class is used

  // constructors
  Block() {}
  Block(const std::string &name, int level) : level_(level), name_(name) {}
  // destructor
  ~Block() { for (auto b:blocks_) delete b; }
  // member functions
  template <typename T> operator std::vector<T>();               // Implicit conversion returning a std::vector of typename T
  std::vector<double> get_row() const { return(values_); }       // Return the values vector
  std::vector<double> get_column(const int n);                   // Return a 1D vector of column n
  std::vector< std::vector<double> > get_2D_matrix();            // Return a 2D vector of column n
  std::vector< std::vector <double> > get_symmetric_2D_matrix(); // Return a symmetric 2D vector of column n
  std::vector<std::string> & get_names(void);                    //
  int nrows(void) const { return int(blocks_.size());}           // Return number of rows
  int ncols(void) const { return int(values_.size());}           // Return number of columns
  int nrows_not_match(const std::string &pattern) { return nrows()-nrows_match(pattern);}
  int nrows_match(const std::string &pattern);
  Block & operator[](const char *key);
  Block & operator[](const std::string &keyword);
  //template <typename T>
  double operator[](const int ind) {
    if (int(values_.size())<ind+1) {
      std::cerr << std::endl << "ERROR reading input-file: '" << name_ << "': index outside limit, " << ind << " > " << values_.size()-1 << std::endl << std::endl;
      std::cerr << values_[0] << std::endl;
      exit(1);
    }
    return values_[ind];
  }
  operator int() const { return int(values_[0]); }
  //operator std::vector<double>(){ return get_vector<double>(); }
  //operator std::vector<int>(){ return get_vector<int>(); }
  operator double(){ return values_[0]; }
  //std::ostream& operator<<(std::ostream& out){ return out << values[0];};
  void print(void) {
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
  void print_value(const std::string &keyword);

private:
    void remove_key_identifyer(std::string &start_id, std::string &end_id);
    Block * find(const std::string &name);
    Block * find_or_create(const std::string &name);
    template <typename T>
    void print_vec(std::vector<T> &vec) {
      for (const auto& v:vec)
        std::cout << v << " ";
      std::cout << std::endl;
    }
    //void print_ivec(std::vector<int> &v);
    //void print_svec(std::vector<std::string> &v);

    friend class Input;
};

//--------------------------------------------------------------
//
//   Declarations of template functions
//
//--------------------------------------------------------------

// Implicit conversion returning a std::vector of typename T
// Example: std::vector<int> size = input["size-vector"];
template <typename T>
Block::operator std::vector<T> () {
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



//template <typename T>
class Input {
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
    // constructors
    Input(const std::string filename)
        : head_block(new Block())
    {
        end_word = key_start_id + end_word + key_end_id;
        current_block.push(head_block);
        set_input_from_file(filename);
        process_input();
    }

    // destructor
    ~Input() {
        current_block.pop();
        delete head_block;
    }

    void print(void) { head_block->print(); }
    Block & operator[](const char *key);


private:
    void process_input();
    void set_input_from_file(const std::string &filename);
    void read_set(std::istringstream *stream);
    void read_keyword(std::string &word, std::istringstream *stream);
    void read_end(void);
    void read_block_content(std::string &word, std::istringstream *stream);
    void read_block_content_v2(std::string &word, std::istringstream *stream);
    void remove_key_identifiers(void) { head_block->remove_key_identifyer(key_start_id, key_end_id); }
    bool is_keyword(const std::string &word);
    bool is_numeric (const std::string& str) const;
    int open(const std::string filename);
    void remove_space(std::string &str);
    void remove_comments(std::string &str);
    const std::string inc_string(std::string &s) const { return std::to_string(std::stoi(s)+1); }
    const std::string& remove_key_tags(std::string& name);
    template <typename T>
    void push_back_word(const std::string& word, Block* newblock) {
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
