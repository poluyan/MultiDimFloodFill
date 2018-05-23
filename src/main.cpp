#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <set>

template <typename T>
T pdf(std::vector<T> x)
{
    T rez = 0;
    for(auto i : x)
        rez += std::pow(i, 2.0);

    return std::sqrt(rez) < 2 ? 1 : 0;
}

void FloodFill_MultipleGrids_VonNeumann(std::vector<std::vector<double>>& grids,
                                        std::vector<std::vector<int>> &points,
                                        std::set<std::vector<int>> &visited,
                                        std::vector<std::vector<double>> &samples,
                                        std::vector<double> dx,
                                        size_t &counter,
                                        size_t &fe_count)
{
    while(!points.empty())
    {
        auto t = points.back();
        points.pop_back();
        
        auto it = visited.find(t);
        if(!(it == visited.end()))
        {
            counter++;
            continue;
        }
        visited.insert(t);

        std::vector<double> dot(t.size());
        for(size_t i = 0; i != dot.size(); i++)
        {
            dot[i] = grids[i][t[i]] + dx[i];
        }

        double val = pdf(dot);
        fe_count++;
        if(val > 0.5)
        {
            std::vector<int> point = t;
            samples.push_back(dot);

            // n-dimensional Manhattan distance with r = 1
            for(size_t i = 0; i != t.size(); i++)
            {
                point = t;
                point[i] = point[i] + 1;

                if(point[i] < 0 || point[i] > grids[i].size() - 1)
                    continue;

                points.push_back(point);
            }
            for(size_t i = 0; i != t.size(); i++)
            {
                point = t;
                point[i] = point[i] - 1;

                if(point[i] < 0 || point[i] > grids[i].size() - 1)
                    continue;

                points.push_back(point);
            }
        }
    }
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
        std::vector<double> onegrid(num_points);
        double startp = lb;
        double endp = ub;
        double es = endp - startp;
        for(size_t j = 0; j != onegrid.size(); j++)
        {
            onegrid[j] = startp + j*es/(num_points);
        }
        grids[i] = onegrid;
        dx[i] = es/(num_points*2);
    }

    // finding start dot over created grid
    std::vector<int> startdot(init_point.size());
    for(size_t i = 0; i != startdot.size(); i++)
    {
        std::vector<double> val(grids[i].size());
        for(size_t j = 0; j != val.size(); j++)
        {
            val[j] = grids[i][j] + dx[i];
        }
        auto pos1 = std::lower_bound(val.begin(), val.end(), init_point[i]);
        startdot[i] = std::distance(val.begin(), pos1) - 1;
    }

    std::vector<std::vector<int>> points;
    std::vector<std::vector<int>> boundary;

    points.push_back(startdot);

    std::set<std::vector<int>> visited;

    std::vector<std::vector<double> > samples;

    size_t counter = 0;
    size_t fe_count = 0;

    FloodFill_MultipleGrids_VonNeumann(grids, points, visited, samples, dx, counter, fe_count);

    std::cout << counter << std::endl;
    std::cout << "fe count: " << fe_count << std::endl;
    std::cout << "samples: " << samples.size() << std::endl;
    std::cout << samples.size()/double(fe_count) << std::endl;

}

int main()
{
    std::vector<double> start = {0, 0};
    b4MultipleGrids(start);
}
