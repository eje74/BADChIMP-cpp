#include <vector>
#include <array>
#include <iostream>

template <typename T>
using vec = std::vector<T>;
//template <typename T>
//using mat = vec<std::vector<T>>;

struct voxel
{
    static constexpr int type = 11;
    static constexpr int dim = 3;
    static constexpr int n = 8;
    static constexpr std::array<std::array<int, dim>, n> vertex = {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}};    
};

template <typename T>
class Cell
{
private:
    std::vector<double> points_;
    std::vector<int> conn_;

public:
    Cell(std::vector<std::vector<double>> &nodes, std::vector<int> &size)
        : points_(), conn_()
    {
        conn_.reserve(T::n * nodes.size());
        std::vector<int> index_list(1 + size[0] * size[1] * size[2], -1);
        for (const auto &n : nodes)
        {
            for (const auto &v : T::vertex)
            {
                std::vector<double> p(T::dim);
                for (auto i = 0; i < T::dim; i++)
                    p[i] = n[i] + v[i];
                // give cell-points a unique index
                int idx = p[0] + p[1] * size[0] + p[2] * size[0] * size[1];
                if (index_list[idx] < 0)
                {
                    points_.insert(points_.end(), std::begin(p), std::end(p));
                    index_list[idx] = int(points_.size() / T::dim) - 1;
                }
                conn_.push_back(index_list[idx]);
            }
        }
    };
};

