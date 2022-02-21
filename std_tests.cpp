#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <execution>

// template <typename T, typename U>
// U dublett(const T&)
template <typename T, typename U>
void zip(const T a_beg, const T a_end, const T b_beg, U ab_beg)
//-----------------------------------------------------------------------------------
{    
    std::cout << sizeof(a_beg) << std::endl;
    std::vector<int> nr(std::distance(a_beg, a_end));
    std::iota(std::begin(nr), std::end(nr), 0);    
    std::for_each(std::execution::par_unseq, std::begin(nr), std::end(nr), [&](const auto& i){
        *std::next(ab_beg, i*2)   = *std::next(a_beg, i);
        *std::next(ab_beg, i*2+1) = *std::next(b_beg, i);
    });
}

template <typename T>
int count(T a_beg, T a_end) {
    return std::distance(a_beg, a_end);
}

int main() 
{
    // std::vector<int> A {0,0,2,2,4,4,6,6};
    // std::vector<int> B {1,1,3,3,5,5,7,7};
    const std::vector<int> A {0,2,4,6};
    const std::vector<int> B {1,3,5,7};
    std::vector<int> AB(2 + std::size(A) + std::size(B), 0);
    zip(std::begin(A), std::end(A), std::begin(B), std::begin(AB)+1);
    
    auto it = std::begin(A);
    std::cout << count(it, std::end(A)) << std::endl;
    // std::vector<int> nr(std::size(AB));
    // std::iota(std::begin(nr), std::end(nr), 0);

    // // std::vector<int> AB(std::size(A));
    // // int i = 0;
    // auto a = std::begin(A);
    // auto b = std::begin(B);
    // auto ab = std::begin(AB);
    // std::for_each(std::execution::par_unseq, std::begin(nr), std::end(nr), [&](const auto& i){
    //     *std::next(ab, i*2) = *std::next(a, i);
    //     *std::next(ab, i*2+1) = *std::next(b, i);
    // });
    // std::transform(std::begin(A), std::end(A), std::begin(B), std::begin(AB), [](int a, int b){
    //     //std::swap(a,b);
    //     return i++%2 ? b : a;
    // });
    //std::copy(std::begin(AB), std::end(AB), std::ostream_iterator<int>(std::cout, ", "));
    for (auto& ab : AB)
        std::cout << ab << ", ";
    std::cout << std::endl;
}
