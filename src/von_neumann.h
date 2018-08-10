#ifndef VON_NEUMANN_H_INCLUDED
#define VON_NEUMANN_H_INCLUDED

#include <vector>
#include <set>
#include <cmath>
#include <iostream>
#include "trie.h"

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

void FloodFill_MultipleGrids_VonNeumann(std::vector<std::vector<double>>& grids,
                                        std::vector<std::vector<int>> &points,
                                        std::set<std::vector<int>> &visited,
                                        std::vector<std::vector<double>> &samples,
                                        std::vector<double> dx,
                                        size_t &counter,
                                        size_t &fe_count);

void FloodFill_MultipleGrids_VonNeumann_trie(std::vector<std::vector<double>>& grids,
        std::vector<int> &start,
        trie_cpp::Trie<trie_cpp::Node<int>,int> &samples,
        std::vector<double> dx,
        size_t &counter,
        size_t &fe_count);

void b4MultipleGrids(std::vector<double> init_point);
void b4MultipleGrids_trie(std::vector<double> init_point);

#endif
