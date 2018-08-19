#include <iostream>
#include <algorithm>
#include <vector>
#include "timer.h"
#include "trie.h"
#include "print2file.h"
#include "von_neumann.h"

int main()
{
    size_t dim = 3;
    size_t grid_sizes = 20;
    
    std::vector<double> start(dim, 0); // std::sqrt(0.5)
    start.front() = -2;

    timer::Timer time_cpp11;

    time_cpp11.reset();
    b4MultipleGrids(start, grid_sizes, true);
    std::cout << "total time:        " << time_cpp11.elapsed_seconds() << std::endl;
    time_cpp11.reset();
    b4MultipleGrids_trie(start, grid_sizes, true);
    std::cout << "total time (trie): " << time_cpp11.elapsed_seconds() << std::endl;
}
