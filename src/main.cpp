#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <set>
#include <fstream>
#include <string>
#include "timer.h"
#include "trie.h"

size_t grid_sizes = 20;
size_t dimension = 2;

int alpha_size = grid_sizes + 1;

template <typename T>
bool pdf(std::vector<T> x)
{
//    T rez1 = 0, rez2 = 0;
//    for(auto i : x)
//    {
//        rez1 += std::pow(i+0.5, 2.0);
//        rez2 += std::pow(i-0.5, 2.0);
//    }
//    return (std::sqrt(rez1) < 1.5 || std::sqrt(rez2) < 1.5) && std::sqrt(rez2) > 0.5 ? true : false;
    
//    T sum1 = 0, sum2 = 0;
//    for(auto i : x)
//    {
//        sum1 += i*i;
//        sum2 += (i+1)*(i+1);
//    }
//
//    return std::exp(-5.0*std::abs(sum1 - 1.0)) > 0.1 && std::sqrt(sum2) > 1.75 ? true : false ;

    T rez = 0;
    if(x.size() == 2)
    {
        T a = 2.0, b = 0.5;
        rez = std::signbit(std::pow(x[0]*x[0]+x[1]*x[1]-(a*a+b*b),2.0)-4.0*a*a*(b*b));
    }
    if(x.size() == 3)
    {
        T a = 2.0, b = 0.5;
        rez = std::signbit(std::pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]-(a*a+b*b),2.0)-4.0*a*a*(b*b-x[2]*x[2]));
    }
    return rez;
}

template <typename T>
void print2file2d(std::string fname, std::vector<std::vector<T> > u)
{
    std::ofstream fOut;
    fOut.open(fname.c_str());
    if(!fOut.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    fOut.precision(10);
    for(size_t i = 0; i < u.size(); i++)
    {
        for(size_t j = 0; j < u[i].size(); j++)
        {
            fOut << std::scientific << u[i][j] << '\t';
        }
        fOut << std::endl;
    }
    fOut.close();
    std::cout << fname << std::endl;
}

template <typename T>
void print2file_3d(std::string fname, std::vector<std::vector<T> > u)
{
    std::ofstream fOut;
    fOut.open(fname.c_str());
    if(!fOut.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }
    fOut.precision(10);
    for(size_t i = 0; i < u.size(); i++)
    {
        for(size_t j = 0; j < u[i].size(); j++)
        {
            fOut << std::scientific << u[i][j] << '\t';
        }
        fOut << 1.0 << std::endl;
    }
    fOut.close();
    std::cout << fname << std::endl;
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

        std::string s(t.begin(), t.end());
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

void FloodFill_MultipleGrids_VonNeumann_trie(std::vector<std::vector<double>>& grids,
        std::vector<std::vector<int>> &points,
        trie::TrieNode *visited,
        std::vector<std::vector<double>> &samples,
        std::vector<double> dx,
        size_t &counter,
        size_t &fe_count)
{
    while(!points.empty())
    {
        auto t = points.back();
        points.pop_back();

        if(trie::search(visited,t))
        {
            counter++;
            continue;
        }
        trie::insert(visited, t, alpha_size);

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

void b4MultipleGrids_trie(std::vector<double> init_point)
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

    trie::TrieNode *visited = trie::getNode(alpha_size);

    std::vector<std::vector<double> > samples;

    size_t counter = 0;
    size_t fe_count = 0;

    FloodFill_MultipleGrids_VonNeumann_trie(grids, points, visited, samples, dx, counter, fe_count);

    std::cout << counter << std::endl;
    std::cout << "fe count: " << fe_count << std::endl;
    std::cout << "samples: " << samples.size() << std::endl;
    std::cout << samples.size()/double(fe_count) << std::endl;

    if(dim == 2) print2file2d("maps/sample2d.dat", samples);
    else if(dim == 3) print2file_3d("maps/sample3d.dat", samples);
}



int main()
{
    size_t dim = dimension;
    std::vector<double> start(dim, 0); // std::sqrt(0.5)
    start.front() = -2;

    Timer time_cpp11;

    time_cpp11.reset();
    b4MultipleGrids(start);
    std::cout << time_cpp11.elapsed_seconds() << std::endl;
    time_cpp11.reset();
    b4MultipleGrids_trie(start);
    std::cout << time_cpp11.elapsed_seconds() << std::endl;
}
