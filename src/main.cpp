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
    size_t dim = 3;
    size_t grid_sizes = 20;
    
    std::vector<double> start(dim, 0); // std::sqrt(0.5)
    start.front() = -2;

    timer::Timer time_cpp11;

    time_cpp11.reset();
    b4MultipleGrids_VonNeumann(start, grid_sizes, true);
    std::cout << "VonNeumann total time:        " << time_cpp11.elapsed_seconds() << std::endl;
    
    
    time_cpp11.reset();
    b4MultipleGrids_VonNeumann_deque(start, grid_sizes, true, 2);
    std::cout << "VonNeumann deque total time:        " << time_cpp11.elapsed_seconds() << std::endl;
        
    time_cpp11.reset();
    b4MultipleGrids_VonNeumann_trie(start, grid_sizes, true, true, 2);
    std::cout << "VonNeumann total time (trie): " << time_cpp11.elapsed_seconds() << std::endl;
    
    time_cpp11.reset();
    b4MultipleGrids_Moore(start, grid_sizes, true);
    std::cout << "Moore total time:        " << time_cpp11.elapsed_seconds() << std::endl;
}
