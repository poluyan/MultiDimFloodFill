#ifndef VON_NEUMANN_H
#define VON_NEUMANN_H

#include <vector>
#include <set>
#include <cmath>
#include <deque>
#include <iostream>
#include "trie.h"

template <typename T, typename V>
struct cell
{
    std::vector<T> dot;
    V value;

    cell() {}
    cell(std::vector<T> _dot, V _value):dot(_dot),value(_value) {}
    cell(const cell &a):dot(a.dot), value(a.value) {}
    bool operator==(const cell& a) const
    {
        return dot == a.dot;
    }
    bool operator<(const cell& a) const
    {
        return dot < a.dot;
    }
};

void FloodFill_MultipleGrids_VonNeumann(
    const std::vector<std::vector<double>> &grids,
    std::vector<std::vector<int>> &points,
    std::set<std::vector<int>> &visited,
    std::vector<std::vector<double>> &samples,
    const std::vector<double> &dx,
    size_t &counter,
    size_t &fe_count,
    bool outside_bounds);

void FloodFill_MultipleGrids_VonNeumann_deque(
    const std::vector<std::vector<double>> &grids,
    std::deque<cell<int, int>> &points,
    std::set<std::vector<int>> &visited,
    std::vector<std::vector<double>> &samples,
    const std::vector<double> &dx,
    size_t &counter,
    size_t &fe_count,
    bool outside_bounds,
    int thread_number);

void FloodFill_MultipleGrids_VonNeumann_trie(
    const std::vector<std::vector<double>> &grids,
    std::vector<int> &start,
    trie_cpp::Trie<trie_cpp::Node<int>,int> &samples,
    const std::vector<double> &dx,
    size_t &counter,
    size_t &fe_count,
    bool outside_bounds,
    bool multithread,
    int thread_number);

void b4MultipleGrids_VonNeumann(const std::vector<double> &init_point,
                                size_t grid_sizes,
                                bool outside_bounds);

void b4MultipleGrids_VonNeumann_deque(const std::vector<double> &init_point,
                                      size_t grid_sizes,
                                      bool outside_bounds,
                                      int thread_number);

void b4MultipleGrids_VonNeumann_trie(const std::vector<double> &init_point,
                                     size_t grid_sizes,
                                     bool outside_bounds,
                                     bool multithread,
                                     int thread_number);

#endif
