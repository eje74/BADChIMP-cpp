/*
 * Input.cpp
 *
 *  Created on: 4. mar. 2015
 *      Author: janlv
 */

#include "Input.h"

//---------------------
// operator []
//---------------------
Block & Block::operator[](const char *key) {
  if (blocks.empty()) {
    return *this;
  } else {
    std::string keyword(key);
    //keyword.assign(key);
    Block *ret = find(keyword);
    if (ret==nullptr) {
      std::cerr << "ERROR! Missing keyword in <"<< name << "> block: \'" << keyword << "\' not found!" << std::endl;
      exit(1);
      // Could use boost.optional here : boost::optional<Block &> Block::operator[], and return boost::optional<Block &>();
    }
    return *ret;
  }
}

//---------------------
// operator []
//---------------------
Block & Block::operator[](const std::string &keyword) {
  if (blocks.empty()) {
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

//---------------------
// operator []
//---------------------
double Block::operator[](const unsigned int ind) {
  if (values.size()<ind+1) {
    std::cerr << std::endl << "ERROR reading input-file: '" << name << "' index " << ind << " is outside array limits " << values.size()-1 << std::endl << std::endl;
    exit(1);
  }
  return values[ind];
};

//--------------------------------------------------
// create and return a vector of column n of a block
//--------------------------------------------------
std::vector<double> Block::get_column(int n) {
  std::vector<double> col;
  for (int i=0; i<(int)blocks.size(); i++) {
    if ( (n+1)>blocks[i]->values.size()) {
      std::cerr << "ERROR in Block::get_column: index " << n
          << " beyond limit " << blocks[i]->values.size()-1 << std::endl;
      exit(1);
    }
    col.push_back(blocks[i]->values[n]);
  }
  return col;
};

//--------------------------------------------------
// create and return a MATRIX of size nrows by ncols
//--------------------------------------------------
std::vector< std::vector <double> > Block::get_matrix() {
  int nrows = this->nrows();
  int ncols = (*this)[std::to_string(0)].ncols();
  std::vector< std::vector<double> > mat(nrows, std::vector<double>(ncols));
  int n = 0;
  for (auto &row : mat) {
    row = (*this)[std::to_string(n++)].get_row();
  }
  return mat;
};

//--------------------------------------------------
// create and return a SYMMETRIC MATRIX of size nrows by ncols
// where the lower triangle is overwritten by the upper triangle
//--------------------------------------------------
std::vector< std::vector <double> > Block::get_symmetric_matrix() {
  std::vector< std::vector<double> > mat = get_matrix();
  int nrows = this->nrows();
  int ncols = (*this)[std::to_string(0)].ncols();
  if (nrows != ncols) {
    std::cerr << "ERROR! The matrix is not quadratic (" << nrows << ", " << ncols << "), unable to create symmetric matrix!" << std::endl;
    exit(1);
  }
  for (int i=0; i<nrows; ++i) {
    for (int j=0; j<ncols; ++j) {
      if (j>i) {
        mat[j][i] = mat[i][j];
      }
    }
  }
  return mat;
}


//---------------------
// print
//---------------------
void Block::print_value(const std::string &keyword) {
  for (int i=0; i<(int)values.size(); i++) {
    std::cout << values[i] << " ";
  }
  std::cout << std::endl;
}

//---------------------
// print
//---------------------
void Block::print(void) {
  std::string indent = std::string(level*3, ' ');
  std::cout << indent+name << ":";
  for (auto &val : values) {
    std::cout << val << " ";
  }
  std::cout << std::endl;

  // recursive call
  for (auto &block : blocks) {
    block->print();
  }

}

//---------------------
// return number of rows matching given string
//---------------------
int Block::nrows_match(const std::string &pattern) {
  int m=0;
  for (int i=0; i<(int)blocks.size(); i++) {
    if (blocks[i]->name == pattern)
      m++;
  }
  // needs re-implementation
  //  for (int i=0; i<(int)var_names.size(); i++) {
  //    if (var_names[i].find(pattern)<std::string::npos)
  //      m++;
  //  }
  return m;
}

//---------------------
//
//---------------------
std::vector<std::string> & Block::get_names(void) {
  std::vector<std::string> *var_names = new std::vector<std::string>();
  for (int i=0; i<(int)blocks.size(); i++) {
    var_names->push_back(blocks[i]->name);
  }
  return (*var_names);
}

//---------------------
//
//---------------------
void Block::remove_key_identifyer(std::string &start_id, std::string &end_id) {
  size_t len = name.length();
  if ( (name.substr(0,1)==start_id) && (name.substr(len-1,1)==end_id) ) {
    name = name.substr(1);
    name = name.substr(0, len-2);
  }
  // recursive call for blocks
  for (int i=0; i<(int)blocks.size(); i++)
    blocks[i]->remove_key_identifyer(start_id, end_id);
}

//---------------------
//
//---------------------
Block * Block::find(const std::string &key) {
  //std::cout << "Trying to find " << key << " in block " << name << ".....";
  std::vector<Block*>::iterator it;
  for(it=blocks.begin(); it!=blocks.end(); it++) {
    if ((*it)->name == key) {
      //std::cout << "success!!" << std::endl;
      return *it;
    }
  }
  return nullptr;
}

//---------------------
//
//---------------------
Block * Block::find_or_create(const std::string &name) {
  Block *b = find(name);
  if (b != NULL) return b;
  // did not exist, create new
  //std::cout << "Creating new block with name " << name << std::endl;
  blocks.push_back(new Block(name, level+1));
  return blocks.back();
}

//---------------------
//
//---------------------
void Block::print_dvec(std::vector<double> &v) {
  for (std::vector<double>::const_iterator it=v.begin(); it!=v.end(); it++)
    std::cout << *it << " ";
  std::cout << std::endl;
}

//---------------------
//
//---------------------
void Block::print_ivec(std::vector<int> &v) {
  for (std::vector<int>::const_iterator it=v.begin(); it!=v.end(); it++)
    std::cout << *it << " ";
  std::cout << std::endl;
}

//---------------------
//
//---------------------
void Block::print_svec(std::vector<std::string> &v) {
  for (std::vector<std::string>::const_iterator it=v.begin(); it!=v.end(); it++)
    std::cout << *it << " ";
  std::cout << std::endl;
}

//////////////////////////////
//
//  I N P U T class
//
//////////////////////////////

//---------------------
//
//---------------------
void Input::read(const std::string filename) {
  init(filename);
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
  remove_key_identifiers();
}

//---------------------
// operator []
//---------------------
Block & Input::operator[](const char *key) {
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

//--------------------------------------------
// fill queue of stringstreams from input-file
//-------------------------------------------
void Input::init(const std::string &filename) {
  //std::cout << "Input::init: using file " << filename << std::endl;
  head_block->name = filename;
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


//---------------------
//
//---------------------
void Input::read_set(std::istringstream *stream) {
  double val;
  std::string var;
  (*stream) >> var >> val;
  Block *block = current_block.top()->find_or_create(var);
  block->values.push_back(val);
}

//---------------------
//
//---------------------
void Input::read_keyword(std::string &word, std::istringstream *stream) {
  Block *newblock = current_block.top()->find_or_create(word);  // find the matching state
  current_block.push(newblock);        // add to stack
  (*stream) >> word;                // try to read the next word
  while (stream->good()) {
    if ( is_numeric(word) ) {    // if number: add to range (int)
      //newblock->range.push_back(std::stoi(word)); // std::stoi is C++11
      int num;
      std::istringstream(word) >> num;
      newblock->range.push_back(num);
    } else {
      std::cout << "WARNING! Input ignores non-numeric words after keywords: "
          << "'" << word << "'" << " given after " << current_block.top()->name <<std::endl;
    }
    (*stream) >> word;
  }
}

//---------------------
// If the first word is a number we assume it is a data-value
// and use the line-number of the block to name the block
//
//---------------------
void Input::read_block_content(std::string &word, std::istringstream *stream) {
  Block *parent = current_block.top();
  
  // create name for the new block
  std::string name;
  if ( is_numeric(word) ) {
    // use line-number as name
    name = "0";
    if (parent->blocks.size()>0)
      name = inc_string(parent->blocks.back()->name);
    //std::cout << name << std::endl;
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
      //block->values.push_back(std::stod(word));  // C++11 version
      double num;
      std::istringstream(word) >> num;
      newblock->values.push_back(num);
      //std::cout << "adding  " << num << std::endl;
    } else {
      // create new blocks for each additional word on a line
      newblock = newblock->find_or_create(word);
    }
    if ( !((*stream) >> word) )
      break;
  }
}

//--------------------------------------------------
// increment a numerical string by 1 and return it
//--------------------------------------------------
std::string Input::inc_string(std::string &s) {
  int i;
  std::istringstream(s) >> i;
  i++;
  std::ostringstream os;
  if (os << i) {
    return os.str();
  } else {
    return s;
  }
}


//---------------------
//
//---------------------
bool Input::is_keyword(const std::string &word) {
  size_t a = word.find_first_of(key_start_id);
  if (a < std::string::npos) {
    size_t b = word.find_first_of(key_end_id);
    if (b < std::string::npos) {
      if (a>0 || b<word.length()-1) {
        std::cerr << "Input error! Keyword identifiers found inside keyword: "+word << std::endl;
        exit(-1);
      }
      return true;
    } else {
      std::cerr << "Input error! Missing " + key_end_id + " in keyword: "+word << std::endl;
      exit(-1);
    }
  }
  return false;
}

//---------------------
//
//---------------------
int Input::open(const std::string filename) {
  infile.open(filename.c_str());
  if (!infile) {
    std::cerr << "Error! Could not open file " + filename << std::endl;
    exit(1);
  }
  return 1;
}

//---------------------
//
//---------------------
void Input::remove_space(std::string &str) {
  // remove leading and trailing space
  size_t a = str.find_first_not_of(' ');
  if (a == std::string::npos) a = 0;
  str = str.substr(a, str.find_last_not_of(' ') + 1 - a);
}

//---------------------
//
//---------------------
void Input::remove_comments(std::string &str) {
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

//---------------------
//
//---------------------
bool Input::is_numeric (const std::string& str) {
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


