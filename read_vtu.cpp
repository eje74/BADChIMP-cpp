#include <vector>
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
#include <algorithm>
#include <numeric>
#include <execution>

//#include "src/io/Input.h"

template <typename T>
std::vector<std::vector<uint8_t>> split(const T& begin, const T& end, const std::vector<uint8_t>& delim)
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

template <typename T>
std::vector<std::vector<uint8_t>> cut(const T& begin, const T& end, uint8_t delim)
{
    std::vector<std::vector<uint8_t>> out(std::distance(begin, end));
    std::transform(begin, end, std::begin(out), [&](const auto& vec)->std::vector<uint8_t>{ 
        return {std::begin(vec), std::find(std::begin(vec), std::end(vec), delim)};
    });
    return out;
}

int main() 
{
    //std::string fname = "ascii.vtu";
    std::string fname = "test.inp";
    // std::basic_ifstream<std::uint8_t> file(fname, std::ios::binary);
    // file.unsetf(std::ios::skipws);
    std::ifstream file(fname, std::ios::binary);
    if (!file) {
        std::cerr << std::endl << "ERROR: Unable to open " << fname << std::endl << std::endl;
        return -1;
    }

    const std::vector<uint8_t> data((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));
    // const std::vector<std::uint8_t> data((std::istreambuf_iterator<std::uint8_t>(file)), (std::istreambuf_iterator<std::uint8_t>()));

    auto lines = split(std::begin(data), std::end(data), {'\n','\r'});
    lines = cut(std::begin(lines), std::end(lines), '#');

    for (auto& line : lines) {
        if (line.empty())
            continue;
        // const auto linewords{ split(std::begin(line), std::end(line), {'\t',' '}) };
        //auto linewords = split(std::begin(line), std::end(line), {'#'});
        // const auto linewords{ split(std::begin(line), std::end(line), {'<','>'}) };
        const auto linewords{ split(std::begin(line), std::end(line), {' '}) };
        if (linewords[0][0] == '<')
            std::cout << std::string(std::begin(linewords[0]), std::end(linewords[0])) << std::endl;
        // std::vector<std::string> words;
        // std::transform(std::begin(linewords), std::end(linewords), std::back_inserter(words), [](const auto& word){
        //     return std::string(std::begin(word), std::end(word));
        // });
        // std::for_each(std::begin(words), std::end(words), [](const auto& w){std::cout << '|' << w << '|' << ", ";});
        // std::cout << std::endl;
    }

    // Input inp(fname, {{"keyword","<>"},{"end","</\\w+>$"},{"comment","?"}}, Input::math_off | Input::var_off | Input::skip_data );
    // std::cout << inp << std::endl;
    // auto piece = inp["VTKFile"]["UnstructuredGrid"]["Piece"];
    // std::vector<double> data = piece["Points"]["points"].read_data(12*3);
    // for (const auto& d:data)
    //     std::cout << d;
    // std::cout << std::endl;
}