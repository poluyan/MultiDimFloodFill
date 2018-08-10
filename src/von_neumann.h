#ifndef VON_NEUMANN_H
#define VON_NEUMANN_H

#include <vector>
#include <set>
#include <cmath>
#include <iostream>
#include "trie.h"

void FloodFill_MultipleGrids_VonNeumann(
    std::vector<std::vector<double>>& grids,
    std::vector<std::vector<int>> &points,
    std::set<std::vector<int>> &visited,
    std::vector<std::vector<double>> &samples,
    std::vector<double> dx,
    size_t &counter,
    size_t &fe_count);

void FloodFill_MultipleGrids_VonNeumann_trie(
    std::vector<std::vector<double>>& grids,
    std::vector<int> &start,
    trie_cpp::Trie<trie_cpp::Node<int>,int> &samples,
    std::vector<double> dx,
    size_t &counter,
    size_t &fe_count);

void b4MultipleGrids(std::vector<double> init_point, size_t grid_sizes);
void b4MultipleGrids_trie(std::vector<double> init_point, size_t grid_sizes);

#endif
