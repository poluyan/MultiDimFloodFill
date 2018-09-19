#include "moore.h"
#include "pdf.h"

void FloodFill_MultipleGrids_Moore(const std::vector<std::vector<double>> &grids,
                                   std::vector<std::vector<int>> &points,
                                   std::set<std::vector<int>> &visited,
                                   std::vector<std::vector<double>> &samples,
                                   const std::vector<double> &dx,
                                   size_t &counter,
                                   size_t &fe_count,
                                   bool outside_bounds)
{

    /// generate all permutations
    std::vector<std::vector<int>> variable_values(grids.size(), std::vector<int>(3));
    for(size_t i = 0; i != variable_values.size(); i++)
    {
        variable_values[i][0] = -1;
        variable_values[i][1] = 0;
        variable_values[i][2] = 1;
    }
    std::vector<std::vector<int>> permut = iterate(variable_values);

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

        bool val = pdf(dot);
        fe_count++;
        if(val)
        {
            std::vector<int> point = t;
            samples.push_back(dot);

            // n-dimensional Moore distance with r = 1
            for(size_t i = 0; i != permut.size(); i++)
            {
                point = t;

                if(outside_bounds)
                {
                    for(size_t j = 0; j != point.size(); j++)
                    {
                        point[j] = point[j] + permut[i][j];
                        if(point[j] < 0)
                        {
                            point[j] = grids[j].size() - 1;
                        }
                        if(point[j] > static_cast<int>(grids[j].size() - 1))
                        {
                            point[j] = 0;
                        }
                    }
                }
                else
                {
                    bool flag = true;
                    for(size_t j = 0; j != point.size(); j++)
                    {
                        point[j] = point[j] + permut[i][j];
                        if(point[j] < 0 || point[j] > static_cast<int>(grids[j].size() - 1))
                        {
                            flag = false;
                            break;
                        }
                    }
                    if(!flag)
                        continue;
                }
                points.push_back(point);
            }
        }
    }
}

void b4MultipleGrids_Moore(const std::vector<double> &init_point, size_t grid_sizes, bool outside_bounds)
{
    size_t dim = init_point.size();

    std::vector<double> grid(dim, grid_sizes);
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

    FloodFill_MultipleGrids_Moore(grids, points, visited, samples, dx, counter, fe_count, outside_bounds);

    std::cout << counter << std::endl;
    std::cout << "fe count: " << fe_count << std::endl;
    std::cout << "samples: " << samples.size() << std::endl;
    std::cout << samples.size()/double(fe_count) << std::endl;

    //print2file2d("maps/sample2d.dat", samples);
}
