#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

template <typename T>
T pdf(std::vector<T> x)
{
    T rez = 0;
    for(auto i : x)
    {
        rez += std::pow(i, 2.0);
    }
    return rez;
}

void b4MultipleGrids(std::vector<double> init_point)
{
    size_t dim = init_point.size();

    std::vector<double> gridN = {100, 100};
    std::vector<std::vector<double>> grids(gridN.size());
    std::vector<double> dx(gridN.size());

    for(size_t i = 0; i != grids.size(); i++)
    {
        size_t grid_N = gridN[i];
        std::vector<double> grid(grid_N + 1);
        double startp = -3;
        double endp = 3;
        double es = endp - startp;
        for(size_t j = 0; j != grid.size(); j++)
        {
            grid[j] = startp + j*es/(grid_N);
        }
        grids[i] = grid;
        dx[i] = es/(grid_N*2);
    }
}

int main()
{
    std::vector<double> start = {0, 0};
}
