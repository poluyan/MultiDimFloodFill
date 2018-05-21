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

    // finding start dot over created grid
    std::vector<int> startdot(dim);
    for(size_t i = 0; i != startdot.size(); i++)
    {
        auto pos1 = std::lower_bound(grid.begin(), grid.end(), init_point[i]);
        startdot[i] = std::distance(grid.begin(), pos1) - 1;
    }

    std::vector<std::vector<int>> points;
    std::vector<std::vector<int>> boundary;

    points.push_back(startdot);

    std::set<std::vector<int>> visited;

    std::vector<std::vector<double> > samples;

    size_t counter = 0;
    size_t fe_count = 0;

}

int main()
{
    std::vector<double> start = {0, 0};

}
