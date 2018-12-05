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

//#define _DEBUG_

class Block {
 public:
  int level;                 // used for indenting output
  std::string name;
  //std::vector<int> range;
  //std::vector<std::string> modifiers;
  std::vector<Block*> blocks;
  //std::vector<double> values;
  std::vector<double> values;
  std::string datatype;

  // constructors
  Block(const std::string &name_="", int level_=0) : level(level_), name(name_) {};
  // destructor
  ~Block() { for (int i=0; i<(int)blocks.size(); i++) delete blocks[i]; };

  template <typename T>
  operator std::vector<T> () {
    if (blocks.size()>0) {
      std::vector<T> vec;
      vec.reserve(nrows()*blocks[0]->ncols());
      for (const auto& bl:blocks) {
        vec.insert(vec.end(), bl->values.begin(), bl->values.end());
      }
      return vec;
    } else {
      return(std::vector<T>(values.begin(), values.end()));
    }
  };
  std::vector<double> get_row() { return(values); };
  std::vector<double> get_column(int n);
  std::vector< std::vector<double> > get_2D_matrix();
  std::vector< std::vector <double> > get_symmetric_2D_matrix();
  std::vector<std::string> & get_names(void);
  int nrows(void) { return (int) blocks.size(); };
  int ncols(void) { return((int)values.size());};
  int nrows_not_match(const std::string &pattern) { return nrows()-nrows_match(pattern); };

  int nrows_match(const std::string &pattern);
  Block & operator[](const char *key);
  Block & operator[](const std::string &keyword);
  double operator[](const unsigned int ind);
  operator int() const { return int(values[0]); }
  //operator std::vector<double>(){ return get_vector<double>(); }
  //operator std::vector<int>(){ return get_vector<int>(); }
  operator double(){ return values[0]; }
  //std::ostream& operator<<(std::ostream& out){ return out << values[0];};
  void print(void) {
    std::string indent = std::string(level*3, ' ');
    std::cout << indent+name << ":";
    for (const auto& val : values) {
      std::cout << val << " ";
    }
    std::cout << std::endl;

    // recursive call
    for (const auto& block : blocks) {
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
 };

  // destructor
  ~Input() {
    current_block.pop();
    delete head_block;
  }

  void print(void) { head_block->print(); };
  Block & operator[](const char *key);


 private:
  void process_input();
  void set_input_from_file(const std::string &filename);
  void read_set(std::istringstream *stream);
  void read_keyword(std::string &word, std::istringstream *stream);
  void read_end(void);
  void read_block_content(std::string &word, std::istringstream *stream);
  void read_block_content_v2(std::string &word, std::istringstream *stream);
  void remove_key_identifiers(void) { head_block->remove_key_identifyer(key_start_id, key_end_id); };
  bool is_keyword(const std::string &word);
  const bool is_numeric (const std::string& str);
  int open(const std::string filename);
  void remove_space(std::string &str);
  void remove_comments(std::string &str);
  const std::string inc_string(std::string &s) const { return std::to_string(std::stoi(s)+1); };
  const std::string& remove_key_tags(std::string& name);
  template <typename T>
  void push_back_word(const std::string& word, Block* newblock) {
    T w;
    std::istringstream iss(word);
    while (iss >> w) {
      if (typeid(w) == typeid(char))
        w -= '0';
      newblock->values.push_back(w);
    }
  }
};

#endif /* SRC_INPUT_H_ */
