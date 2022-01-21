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
#include <valarray>
#include <map>


//------------------------------------------------------------
// Demonstration of the Input-class
//
// Contents of file example.file
//
//
// How to access contents of example.file in code:
//
// Input input("example.file");
// double version = input["version"];                       // version = 1.1
// double pi = input["math"]["pi"];                         // pi = 3.12
// std::vector<int> fibo = input["math"]["series"]["fibo"]; // fibo = 2, 3, 5, 7, 11, 13, 17, 19
// std::vector<int> binary = input["binary"];               // binary = 0, 1, 0, 1, 1, 1, 1, 0
//



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

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block() : name_(), blocks_(), values_(), strings_(), datatype_() { }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block(const std::string &name, Block& not_found, bool missing_ok=false) 
        : name_(name), blocks_(), values_(), strings_(), datatype_(), not_found_(&not_found), missing_ok_(missing_ok) { }
    //-----------------------------------------------------------------------------------

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block(const std::string &name, int level, Block& parent) 
        : level_(level), name_(name), blocks_(), values_(), strings_(), datatype_(), parent_(&parent), 
          not_found_(parent.not_found_), missing_ok_(parent.missing_ok_)
    //-----------------------------------------------------------------------------------
    {
        if ( str_func::is_numeric(name_) ) {
            // use line-number as name
            name_ = "0";
            if (parent_)
                name_ = std::to_string(parent_->blocks_.size());
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
        std::cout << ", strings = " << strings_.size() << ", datatype = " << datatype_ << ", parent = " << ((parent_)? ("None") : parent_->name_);
        std::cout << std::endl;                
    }

    //                                     Block
    //-----------------------------------------------------------------------------------
    Block* parent() const { if (parent_) { return (parent_); } else { error("Block " + name_ + " has no parent!"); return nullptr;} }
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
        //std::cout << name_ << std::endl;
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
    // Implicit conversion returning a std::valarray of typename T
    // Example: std::valarray<int> size = input["size-vector"];
    //
    template <typename T> 
    operator std::slice_array<T>() const           
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
    const Block& operator[](const char *keyword) const
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
    void error(const std::string& msg) const
    //-----------------------------------------------------------------------------------
    {
        std::cerr << std::endl << "*** ERROR during reading of input-file ***" << std::endl;
        std::cerr << "Block " << name_ << ": " << msg << std::endl << std::endl;
        exit(1);
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
    std::map<std::string, std::string> map_;

public:
    //                                     Tag 
    //-----------------------------------------------------------------------------------
    Tag() : map_({{"keyword","<>"}, {"end","<end>"}}) { init(); }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    Tag(const std::initializer_list<std::pair<const std::string, std::string>>& map) : map_(map) { init(); }     
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    void init()     
    //-----------------------------------------------------------------------------------
    { 
        // Add defaults for missing tags 
        map_.insert({"keyword",""});
        map_.insert({"end",""});
        map_.insert({"comment","#"});
        map_.insert({"variable","$"});
        if (not key())
            std::cout << std::endl << "!! WARNING No input keyword identifier is provided !!" << std::endl << std::endl; 
    }

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    char key() const { return (map_.at("keyword").empty()) ? 0 : map_.at("keyword")[0]; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    char key_end() const { return (map_.at("keyword").length()==2) ? map_.at("keyword")[1] : 0; }    
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    const std::string& end() const { return map_.at("end"); }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    const std::string& comment() const { return map_.at("comment"); }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    const std::string& var() const { return map_.at("variable"); }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    void var_off() { map_["variable"] = ""; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    const std::string& operator[](const std::string& var) { return map_[var]; }
    //-----------------------------------------------------------------------------------

    //                                     Tag 
    //-----------------------------------------------------------------------------------
    std::string all() const { return comment()+key()+key_end()+var()+end(); }
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

public:
    static constexpr uint8_t default_flag = 0b0000'0000; 
    static constexpr uint8_t math_off     = 0b0000'0001; 
    static constexpr uint8_t missing_ok   = 0b0000'0010;
    static constexpr uint8_t var_off      = 0b0000'0100;
    //static constexpr uint8_t ignore_case  = 0b0000'1000;

    //                                     Input
    //-----------------------------------------------------------------------------------
    Input(const std::string& filename, const Tag& tag=Tag(), const uint8_t flags=default_flag)
        : tag_(tag), not_found_(), head_block_(filename, not_found_, missing_ok&flags), current_block_(&head_block_), math_symbols_() 
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
    //
    void process_input() 
    //-----------------------------------------------------------------------------------
    {
        std::string word;
        std::ifstream infile;
        open(filename(), infile);
        std::string line;
        while ( std::getline(infile, line) ) {
            remove_space(line);
            remove_comments(line);
            if (line.empty()) {
                line_num_++;
                continue;
            }
            //std::cout << "line: |" << line << "|" << std::endl;            
            std::istringstream iss_line(line);
            if ( iss_line >> word ) {                
                //std::cout << *this << std::endl;
                line_num_++;              
                if (word == tag_.end()) {
                    //std::cout << "END: " << word << std::endl;
                    current_block_ = current_block_->parent_;
                } else if (is_keyword(word)) {
                    if (not recursive_ and current_block_->parent_) {
                        current_block_ = current_block_->parent_;
                    }
                    //std::cout << "KEY: " << word << std::endl;
                    read_keyword(word, iss_line);
                } else {
                    // we are inside a block or sub-block
                    //std::cout << "CONTENT: " << word << ", " << iss_line.str() << std::endl;
                    read_block_content(word, iss_line);
                }
            }
        }
        infile.close();
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
        Block& newblock = current_block_->find_or_create(remove_key_tags(word));  // find the matching state
        while (stream >> word) {
            if (not tag_.key_end()) {
                // No end-tag, add word to the block
                read_data(word, newblock);
            } else {
                // Check if we are inside or outside the keyword
                auto pos = word.find(tag_.key_end());
                if (pos != std::string::npos) {
                    // Inside keyword: add word as datatype format
                    newblock.datatype_ = remove_key_tags(word);                
                } else {
                    // Outside keyword: add word as block values
                    read_data(word, newblock);
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
        if (not tag_.var().empty())
            replace_variables_with_values(word);
        if ( str_func::is_numeric(word) ) {
            // Numeric value, add it
            auto dtype = newblock.parent()->datatype_;
            if (dtype == "char") {
                push_back_word<char>(word, newblock);
            } else if (dtype == "int") {
                push_back_word<int>(word, newblock);
            } else {
                push_back_word<double>(word, newblock);
            }
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
                error("Variable " + tag_.var() + var + " is not defined");
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
    void open(const std::string& filename, std::ifstream& file)
    {
        file.open(filename.c_str());
        if (not file) {
            std::cerr << "Error! Could not open file " + filename << std::endl;
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
        // only keep line up to '#' character
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
    //-----------------------------------------------------------------------------------
    const std::string& remove_key_tags(std::string& name)
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

    //                                     Input
    //-----------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------
    template <typename T>
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
