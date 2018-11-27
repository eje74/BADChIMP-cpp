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

//#define _DEBUG_

class Block {
 public:
  int level;                 // used for indenting output
  std::string name;
  std::vector<int> range;
  std::vector<Block*> blocks;
  std::vector<double> values;

  // constructors
  Block(const std::string &name_="", int level_=0) : level(level_), name(name_) {};
  // destructor
  ~Block() { for (int i=0; i<(int)blocks.size(); i++) delete blocks[i]; };

  std::vector<double> vector() { return(values); };
  std::vector<double> get_row() { return(values); };
  std::vector<double> get_column(int n);
  std::vector< std::vector<double> > get_matrix();
  std::vector< std::vector <double> > get_symmetric_matrix();
  std::vector<std::string> & get_names(void);
  int nrows(void) { return (int) blocks.size(); };
  int ncols(void) { return((int)values.size());};
  int nrows_not_match(const std::string &pattern) { return nrows()-nrows_match(pattern); };

  int nrows_match(const std::string &pattern);
  Block & operator[](const char *key);
  Block & operator[](const std::string &keyword);
  double operator[](const unsigned int ind);
  operator int() const { return int(values[0]); }
  operator std::vector<double>() const { return(values); }
  operator std::vector<int>() const { return(std::vector<int>(values.begin(), values.end())); }
  operator double(){ return values[0]; }
  //std::string get_string(){ return name; }
  //operator int(){ if (level>-1) return values[0]; else {std::cout << name << std::endl; exit(1);} };
  //operator double(){ if (level>-1) return values[0]; else {std::cout << name << std::endl; exit(1);} };
  operator std::vector<double>&(){ return values; };
  //std::ostream& operator<<(std::ostream& out){ return out << values[0];};
  void print(void);
  void print_value(const std::string &keyword);

 private:
  void remove_key_identifyer(std::string &start_id, std::string &end_id);
  Block * find(const std::string &name);
  Block * find_or_create(const std::string &name);
  void print_dvec(std::vector<double> &v);
  void print_ivec(std::vector<int> &v);
  void print_svec(std::vector<std::string> &v);

  friend class Input;
};


class Input {
 private:
  std::string key_start_id, key_end_id, end_word, set_word;
  std::ifstream infile;
  std::stack<Block*> current_block;
  Block *head_block;
  std::queue<std::istringstream*> input;

 public:
  // constructors
  Input(const std::string &start_id="<", const std::string &end_id=">" ,
      const std::string &end="end", const std::string &set="set")
 : key_start_id(start_id),
   key_end_id(end_id),
   end_word(end),
   set_word(set),
   head_block(new Block())
 {
    end_word = key_start_id + end_word + key_end_id;
    current_block.push(head_block);
 };

  // destructor
  ~Input() {
    current_block.pop();
    delete head_block;
  }

  void read(const std::string filename);
  void print(void) { head_block->print(); };
  Block & operator[](const char *key);

 private:
  void init(const std::string &filename);
  void read_set(std::istringstream *stream);
  void read_keyword(std::string &word, std::istringstream *stream);
  void read_end(void);
  void read_block_content(std::string &word, std::istringstream *stream);
  void read_block_content_v2(std::string &word, std::istringstream *stream);
  void remove_key_identifiers(void) { head_block->remove_key_identifyer(key_start_id, key_end_id); };
  bool is_keyword(const std::string &word);
  int open(const std::string filename);
  void remove_space(std::string &str);
  void remove_comments(std::string &str);
  bool is_numeric (const std::string& str);
  std::string inc_string(std::string &s);
};

#endif /* SRC_INPUT_H_ */
