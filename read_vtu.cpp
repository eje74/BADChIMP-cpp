#include <vector>
#include <algorithm>
#include <string_view>
#include <numeric>
#include <execution>
//#include <queue>

#include "src/io/Input.h"

template <typename T>
std::vector<std::vector<uint8_t>> split(const T& begin, const T& end, const std::vector<uint8_t>& delim)
{
    const auto par{ std::execution::par_unseq };

    const std::vector<uint8_t> data(begin, end);
    std::vector<int> nr(std::size(data));
    std::iota(std::begin(nr), std::end(nr), 0);

    std::vector<int> delim_ind(std::size(data));
    std::transform(par, std::begin(data), std::end(data), std::begin(nr), std::begin(delim_ind), [&](const auto& vec, const auto& ind){
        const bool found {std::find(std::begin(delim), std::end(delim), vec) == std::end(delim)};
        return found ? -1 : ind;   
    });

    std::vector<int> word_end_ind;
    std::copy_if(std::begin(delim_ind), std::end(delim_ind), std::back_inserter(word_end_ind), [](const auto& ind){ return ind >= 0; });
    word_end_ind.push_back(std::size(data));

    std::vector<int> word_len_tmp;
    std::adjacent_difference(std::begin(word_end_ind), std::end(word_end_ind), std::back_inserter(word_len_tmp));
    std::vector<int> word_len(word_len_tmp);
    std::transform(std::begin(word_len_tmp)+1, std::end(word_len_tmp), std::begin(word_len)+1, [](const auto& len){ return len-1; });

    std::vector<std::vector<uint8_t>> words_tmp(std::size(word_len));
    std::transform(par, std::begin(word_end_ind), std::end(word_end_ind), std::begin(word_len), std::begin(words_tmp), [&data](const auto& ind, const auto& len){
        return std::vector<uint8_t>(std::begin(data)+ind-len, std::begin(data)+ind);
    }); 

    std::vector<std::vector<uint8_t>> words;
    std::copy_if(std::begin(words_tmp), std::end(words_tmp), std::back_inserter(words), [](const auto& word){ return std::size(word) > 0; });

    return words;
}

int main() 
{
    std::string fname = "0000_fluid_0000000.vtu";
    std::ifstream file(fname, std::ios::binary);
    if (!file) {
        std::cerr << std::endl << "ERROR: Unable to open " << fname << std::endl << std::endl;
        return -1;
    }

    const std::vector<uint8_t> file_data((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));

    const auto lines{ split(std::begin(file_data), std::end(file_data), {'\n','\r'}) };
    for (auto& line : lines) {
        const auto words{ split(std::begin(line), std::end(line), {'\t',' '}) };
        auto word = words.begin();
        while (word != std::end(words)) {
            std::cout << std::string(std::begin(*word), std::end(*word)) << ", ";
            word = std::next(word);
        }
        std::cout << std::endl;
    }

    // Input inp(fname, {{"keyword","<>"},{"end","</\\w+>$"},{"comment","?"}}, Input::math_off | Input::var_off | Input::skip_data );
    // std::cout << inp << std::endl;
    // auto piece = inp["VTKFile"]["UnstructuredGrid"]["Piece"];
    // std::vector<double> data = piece["Points"]["points"].read_data(12*3);
    // for (const auto& d:data)
    //     std::cout << d;
    // std::cout << std::endl;
}