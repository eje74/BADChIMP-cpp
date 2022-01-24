#include <variant>
#include <string>
#include <iostream>
#include <vector>


template <typename T>
class Block
{
    public:
    std::vector<T> data_;
    //std::vector<Block<T>> blocks_;
    Block(const T& val) : data_(val) {}
    friend std::ostream& operator<<(std::ostream& out, const Block& block) { out << "Size: " << block.data_.size() << std::endl; return out; }

};

// struct My_visitor {
//     int operator()(int i) { return i; }
//     double operator()(double d) { return d; }
//     std::string& operator()(const std::string& s) { return s; }
//     // void operator()(int i) { std::cout << "int, " << i << std::endl; }
//     // void operator()(double d) { std::cout << "double, = " << d << std::endl; }
//     // void operator()(const std::string& s) { std::cout << "string, = " << s << std::endl; }
// };
// struct My_visitor {
//     void operator()(int i) { std::cout << "int, " << i << std::endl; }
//     void operator()(const Block<double>& d) { std::cout << "double, = " << d << std::endl; }
//     void operator()(const Block<char>& s) { std::cout << "char, = " << s << std::endl; }
// };

template<class T>
struct streamer {
    const T& val;
};
template<class T> streamer(T) -> streamer<T>;

template<class T>
std::ostream& operator<<(std::ostream& os, streamer<T> s) {
    os << s.val;
    return os;
}

template<class... Ts>
std::ostream& operator<<(std::ostream& os, streamer<std::variant<Ts...>> sv) {
   std::visit([&os](const auto& v) { os << streamer{v}; }, sv.val);
   return os;
}

int main() 
{
    std::vector<std::variant<int, double, std::string>> var_vec;
    var_vec.emplace_back(1);
    var_vec.emplace_back("text");
    var_vec.emplace_back(10.23);
    for (const auto& v:var_vec)
        std::cout << streamer{v} << " ";
    //     //std::visit(My_visitor{}, v);
    std::cout << std::endl;

}