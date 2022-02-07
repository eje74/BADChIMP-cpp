//#include <vector>
#include <regex>
#include "src/io/Input.h"

// void remove_space(std::string &str) 
// {
//     // remove leading and trailing space
//     size_t a = str.find_first_not_of(' ');
//     if (a == std::string::npos) a = 0;
//     str = str.substr(a, str.find_last_not_of(' ') + 1 - a);
// }

// void remove_comments(std::string &str) 
// {
//     // only keep line up to '#' character
//     size_t a = str.find_first_of("#");
//     if (a < std::string::npos)
//         str = str.substr(0,a);
// }

// void open(const std::string& filename, std::ifstream& file)
// {
//     file.open(filename.c_str());
//     if (not file) {
//         std::cerr << "Error! Could not open file " + filename << std::endl;
//         exit(1);
//     }
// }

// void process(const std::string& fname) 
// {
//     std::string word;
//     std::ifstream infile;
//     int line_num = 0;
//     open(fname, infile);
//     std::string line;
//     while ( std::getline(infile, line) ) {
//         remove_space(line);
//         remove_comments(line);
//         if (line.empty()) {
//             line_num++;
//             continue;
//         }
//         //std::cout << "line: |" << line << "|" << std::endl;            
//         std::istringstream iss_line(line);
//         if ( iss_line >> word ) {                
//             //std::cout << *this << std::endl;
//             line_num++;
//             std::cout << "|" << word << "|" << std::endl;               
//         }
//     }
//     infile.close(); 
// }


int main() 
{
//    std::string end = "</CellData>";
//    std::regex regex("</\\w+>$");
//    std::cout << std::regex_search(end, regex) << std::endl; 
    std::string fname = "output/out1/vtu/0000_fluid_0000001.vtu";
    Input inp(fname, {{"keyword","<>"},{"end","</\\w+>$"},{"comment","?"}}, Input::math_off);
    std::cout << inp << std::endl;

}