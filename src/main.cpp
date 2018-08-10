#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <set>
#include "timer.h"
#include "trie.h"
#include "print2file.h"
#include "von_neumann.h"

void f() {}

int main()
{
    
    size_t dim = 3;
    std::vector<double> start(dim, 0); // std::sqrt(0.5)
    start.front() = -2;

    timer::Timer time_cpp11;

    time_cpp11.reset();
    b4MultipleGrids(start);
    std::cout << time_cpp11.elapsed_seconds() << std::endl;
    time_cpp11.reset();
    b4MultipleGrids_trie(start);
    std::cout << time_cpp11.elapsed_seconds() << std::endl;
}
