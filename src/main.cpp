#include <iostream>
#include <algorithm>
#include <vector>
#include "timer.h"
#include "trie.h"
#include "print2file.h"
#include "von_neumann.h"
#include "moore.h"

int main()
{
    
    /// generate all permutations
    std::vector<std::vector<int>> variable_values(2, std::vector<int>(3));
    for(size_t i = 0; i != variable_values.size(); i++)
    {
        variable_values[i][0] = -1;
        variable_values[i][1] = 0;
        variable_values[i][2] = 1;
    }
    std::vector<std::vector<int>> permut = iterate(variable_values);
    
    for(const auto &i : permut)
    {
        for(const auto &j : i)
        {
            std::cout << j << '\t';
        }
        std::cout << std::endl;
    }

    
//    size_t dim = 3;
//    size_t grid_sizes = 20;
//    
//    std::vector<double> start(dim, 0); // std::sqrt(0.5)
//    start.front() = -2;
//
//    timer::Timer time_cpp11;
//
//    time_cpp11.reset();
//    b4MultipleGrids(start, grid_sizes, true);
//    std::cout << "total time:        " << time_cpp11.elapsed_seconds() << std::endl;
//    time_cpp11.reset();
//    b4MultipleGrids_trie(start, grid_sizes, true);
//    std::cout << "total time (trie): " << time_cpp11.elapsed_seconds() << std::endl;
}
