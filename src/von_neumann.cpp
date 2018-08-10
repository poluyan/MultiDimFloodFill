#include "von_neumann.h"
#include "print2file.h"
#include "pdf.h"

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

        bool val = pdf(dot);
        fe_count++;
        if(val)
        {
            std::vector<int> point = t;
            samples.push_back(dot);

            // n-dimensional Manhattan distance with r = 1
            for(size_t i = 0; i != t.size(); i++)
            {
                point = t;
                point[i] = point[i] + 1;

                if(point[i] < 0 || point[i] > static_cast<int>(grids[i].size() - 1))
                    continue;

                points.push_back(point);
            }
            for(size_t i = 0; i != t.size(); i++)
            {
                point = t;
                point[i] = point[i] - 1;

                if(point[i] < 0 || point[i] > static_cast<int>(grids[i].size() - 1))
                    continue;

                points.push_back(point);
            }
        }
    }
}

void FloodFill_MultipleGrids_VonNeumann_trie(std::vector<std::vector<double>>& grids,
        std::vector<int> &start,
        trie_cpp::Trie<trie_cpp::Node<int>,int> &samples,
        std::vector<double> dx,
        size_t &counter,
        size_t &fe_count)
{
    std::vector<double> dot(grids.size());

    trie_cpp::Trie<trie_cpp::Node<int>,int> visited;
    trie_cpp::Trie<trie_cpp::Node<int>,int> points;
    trie_cpp::Trie<trie_cpp::Node<int>,int> not_coumputed;
    points.insert(start);
    ++fe_count;

    while(!points.empty() || !not_coumputed.empty())
    {
        while(!points.empty())
        {
            auto point = points.get_and_remove_last();
            if(visited.search(point) || samples.search(point))
            {
                continue;
            }
            samples.insert(point);

            auto init_point = point;
            for(size_t i = 0; i != point.size(); i++)
            {
                point = init_point;
                point[i] = point[i] + 1;

                if(point[i] < 0 || point[i] > static_cast<int>(grids[i].size() - 1))
                    continue;

                if(!visited.search(point) && !samples.search(point))
                {
                    not_coumputed.insert(point);
                }
            }
            for(size_t i = 0; i != point.size(); i++)
            {
                point = init_point;
                point[i] = point[i] - 1;

                if(point[i] < 0 || point[i] > static_cast<int>(grids[i].size() - 1))
                    continue;

                if(!visited.search(point) && !samples.search(point))
                {
                    not_coumputed.insert(point);
                }
            }
        }

        size_t number_to_points = 0;
        while(!not_coumputed.empty())
        {
            auto point = not_coumputed.get_and_remove_last();
            std::vector<double> values(point.size());
            for(size_t j = 0; j != values.size(); j++)
            {
                values[j] = grids[j][point[j]] + dx[j];
            }
            bool value = pdf(values);
            ++fe_count;

            if(value)
            {
                ++number_to_points;
                points.insert(point);
            }
            else
            {
                visited.insert(point);
            }
        }
        //std::cout << number_to_points << std::endl;
    }
    counter = visited.get_total_count() + samples.get_total_count();
    visited.remove_tree();
    points.remove_tree();
    not_coumputed.remove_tree();
}

void b4MultipleGrids(std::vector<double> init_point, size_t grid_sizes)
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

    FloodFill_MultipleGrids_VonNeumann(grids, points, visited, samples, dx, counter, fe_count);

    std::cout << counter << std::endl;
    std::cout << "fe count: " << fe_count << std::endl;
    std::cout << "samples: " << samples.size() << std::endl;
    std::cout << samples.size()/double(fe_count) << std::endl;

    //print2file2d("maps/sample2d.dat", samples);
}

void b4MultipleGrids_trie(std::vector<double> init_point, size_t grid_sizes)
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

    //trie_c::TrieNode *visited = trie_c::getNode(alpha_size);

    std::vector<std::vector<double> > samples;
    trie_cpp::Trie<trie_cpp::Node<int>,int> trie_samples;

    size_t counter = 0;
    size_t fe_count = 0;

    FloodFill_MultipleGrids_VonNeumann_trie(grids, startdot, trie_samples, dx, counter, fe_count);

    while(!trie_samples.empty())
    {
        auto point = trie_samples.get_and_remove_last();
        std::vector<double> values(point.size());
        for(size_t j = 0; j != values.size(); j++)
        {
            values[j] = grids[j][point[j]] + dx[j];
        }
        samples.push_back(values);
    }


    std::cout << counter << std::endl;
    std::cout << "fe count: " << fe_count << std::endl;
    std::cout << "samples: " << samples.size() << std::endl;
    std::cout << samples.size()/double(fe_count) << std::endl;

    if(dim == 2) print2file2d("maps/sample2d.dat", samples);
    else if(dim == 3) print2file_3d("maps/sample3d.dat", samples);
}
