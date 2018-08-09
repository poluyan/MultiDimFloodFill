#ifndef PRINT2FILE_H_INCLUDED
#define PRINT2FILE_H_INCLUDED

#include <string>
#include <vector>
#include <fstream>

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
        fOut << std::endl;
    }
    fOut.close();
    std::cout << fname << std::endl;
}

#endif