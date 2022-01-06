#include <vector>
#include <valarray>
#include <iostream>
#include <memory>
//#include "src/io/vector_func.h"
//#include "src/io/VTK.h"

template <typename T>
class data_wrapper
{
    public:
    data_wrapper() { std::cout << "data_wrapper constructor" << std::endl; }
    virtual ~data_wrapper() { std::cout << "data_wrapper destructor" << std::endl; };
    virtual const T at(const int pos) const = 0;
    virtual const T* begin() const = 0;
    virtual const T* end() const = 0;
    virtual const size_t size() const = 0;
  };

template <typename T>
class vec_wrapper : public data_wrapper<T>
{
    private:
        const std::vector<T>& data_;
    public:
        vec_wrapper(const std::vector<T>& data) : data_(data) { std::cout << "vec_wrapper constructor" << std::endl; }
        const T at(const int pos) const { return data_[pos]; }
        const T* begin() const { return &data_[0]; }
        const T* end() const { return &data_[data_.size()]; }
        const size_t size() const { return data_.size(); }

};

template <typename T>
class arr_wrapper : public data_wrapper<T>
{
    private:
        const std::valarray<T>& data_;
    public:
        arr_wrapper(const std::valarray<T>& data) : data_(data) {  std::cout << "arr_wrapper constructor" << std::endl;}
        const T at(const int pos) const { return data_[pos]; }
        // const T* begin() const { return data_.begin(); }
        // const T* end() const { return data_.end(); }
        const T* begin() const { return &data_[0]; }
        const T* end() const { return &data_[data_.size()]; }
        const size_t size() const { return data_.size(); }
};


template <typename T>
void info(const data_wrapper<T>& data)
{
    std::cout << data.at(0) << std::endl;
}

class Variable
{
    private:
        std::vector<int> points_;
    public:
        Variable() : points_({1,2,3,4,5}) { }
        //const std::vector<int>& points() const { return points_; } 
        const vec_wrapper<int> points() const { return vec_wrapper<int>(points_); }
};

template <typename T>
void print_vec(const std::vector<T>& vec)
{
    for (const double v : vec)
        std::cout << v << " ";
    std::cout << std::endl;
}

int main() 
{
    std::vector<double> vec = {1,2,3,4,5,6,7,8,9,10};
    vec_wrapper<double> vdata(vec);
    std::valarray<double> arr = {10,20,30,40,50,60,70,80,90,100};
    arr_wrapper<double> adata(arr);
    const data_wrapper<double>& d1 = vdata;
    const data_wrapper<double>& d2 = adata;
    std::vector<double> cpy(d1.size());
    print_vec(cpy);
    std::copy(d1.begin(), d1.end(), cpy.begin());
    print_vec(cpy);
    std::copy(d2.begin(), d2.end(), cpy.begin());
    print_vec(cpy);
    std::unique_ptr<data_wrapper<double>> d3(new vec_wrapper<double>(vec));
    std::cout << d3->at(3) << std::endl;
}

