#ifndef VON_NEUMANN_H
#define VON_NEUMANN_H

#include <vector>
#include <set>
#include <cmath>
#include <iostream>
#include "trie.h"

void FloodFill_MultipleGrids_VonNeumann(
    const std::vector<std::vector<double>> &grids,
    std::vector<std::vector<int>> &points,
    std::set<std::vector<int>> &visited,
    std::vector<std::vector<double>> &samples,
    const std::vector<double> &dx,
    size_t &counter,
    size_t &fe_count,
    bool outside_bounds);

void FloodFill_MultipleGrids_VonNeumann_trie(
    const std::vector<std::vector<double>> &grids,
    std::vector<int> &start,
    trie_cpp::Trie<trie_cpp::Node<int>,int> &samples,
    const std::vector<double> &dx,
    size_t &counter,
    size_t &fe_count,
    bool outside_bounds);

void b4MultipleGrids_VonNeumann(const std::vector<double> &init_point, size_t grid_sizes, bool outside_bounds);
void b4MultipleGrids_VonNeumann_trie(const std::vector<double> &init_point, size_t grid_sizes, bool outside_bounds);

#endif
