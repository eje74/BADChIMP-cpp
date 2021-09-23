#include <vector>
#include <array>
#include <iostream>

template <typename T>
using vec = std::vector<T>;
template <typename T>
using mat = vec<std::vector<T>>;

struct voxel
{
    static constexpr int type = 11;
    static constexpr int dim = 3;
    static constexpr int n = 8;
    static constexpr std::array<std::array<int, dim>, n> vertex = {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}};
    mat<double> points(mat<double>& nodes, vec<int>& size)
    {
        mat<double> points;
        points.reserve(nodes.size()*n);
        for (const auto &n : nodes)
        {
            for (const auto &v : vertex)
            {
                // give cell-points a unique index
                int idx = (n[0]+v[0]) + (n[1]+v[1])*size[0] + (n[2]+v[2])*size[0]*size[1];

            }
        }
    }
};

int main()
{
    int idx = p[0] + p[1] * dim_[0] + p[2] * dim_[0] * dim_[1];
    std::cout << voxel::type << std::endl;
}