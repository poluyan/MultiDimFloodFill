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

    std::vector<double> grid = {100, 100};
    std::vector<std::vector<double>> grids(grid.size());
    std::vector<double> dx(grid.size());
    
    double lb = -3, ub = 3;

    for(size_t i = 0; i != grids.size(); i++)
    {
        size_t num_points = grid[i];
        std::vector<double> onegrid(num_points + 1);
        double startp = lb;
        double endp = ub;
        double es = endp - startp;
        for(size_t j = 0; j != onegrid.size(); j++)
        {
            onegrid[j] = startp + j*es/(num_points);
        }
        grids[i] = grid;
        dx[i] = es/(num_points*2);
    }
}

int main()
{
    std::vector<double> start = {0, 0};
}
