#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <set>
#include <fstream>
#include <string>
#include "timer.h"

size_t grid_sizes = 11;
size_t dimension = 8;

const int alpha_size = 7;
struct TrieNode
{
    struct TrieNode *children[alpha_size];
    bool isEndOfWord;
};

TrieNode *getNode()
{
    TrieNode *pNode =  new TrieNode;
    pNode->isEndOfWord = false;
    for(int i = 0; i < alpha_size; i++)
        pNode->children[i] = 0;
    return pNode;
}

void insert(struct TrieNode *root, std::vector<int> key)
{
    struct TrieNode *pCrawl = root;
    for(int i = 0; i < key.size(); i++)
    {
        int index = key[i];
        if(!pCrawl->children[index])
            pCrawl->children[index] = getNode();
        pCrawl = pCrawl->children[index];
    }
    pCrawl->isEndOfWord = true;
}

bool search(struct TrieNode *root, std::vector<int> key)
{
    struct TrieNode *pCrawl = root;
    for(int i = 0; i < key.size(); i++)
    {
        int index = key[i];
        if(!pCrawl->children[index])
            return false;
        pCrawl = pCrawl->children[index];
    }
    return (pCrawl != NULL && pCrawl->isEndOfWord);
}

// Ershov habr.com/post/151421/

struct node // структура для представления узлов дерева
{
    int* key;
    int len;
    node* link;
    node* next;
    node(int* x, int n) : len(n), link(0), next(0)
    {
        key = new int [n];
        std::copy(key, key + n, x);
    }
    ~node()
    {
        delete[] key;
    }
};

int prefix(int* x, int n, int* key, int m) // длина наибольшего общего префикса строк x и key
{
    for(int k=0; k<n; k++)
        if(k==m || x[k]!=key[k])
            return k;
    return n;
}

node* find(node* t, int* x, int n=0) // поиск ключа x в дереве t
{
    if(!n) n = sizeof(x)/sizeof(x[0])+1;
    if(!t) return 0;
    int k = prefix(x,n,t->key,t->len);
    if(k==0) return find(t->next,x,n);   // поищем у сестры
    if(k==n) return t;
    if(k==t->len) return find(t->link,x+k,n-k);   // поищем у старшей дочери
    return 0;
}

void split(node* t, int k) // разбиение узла t по k-му символу ключа
{
    node* p = new node(t->key+k,t->len-k);
    p->link = t->link;
    t->link = p;
    int* a = new int[k];
    std::copy(a,a + k,t->key);
    delete[] t->key;
    t->key = a;
    t->len = k;
}

node* insert(node* t, int* x, int n=0) // вставка ключа x в дерево t
{
    if(!n) n = sizeof(x)/sizeof(x[0])+1;
    if(!t) return new node(x,n);
    int k = prefix(x,n,t->key,t->len);
    if(k==0) t->next = insert(t->next,x,n);
    else if(k<n)
    {
        if(k<t->len)   // режем или нет?
            split(t,k);
        t->link = insert(t->link,x+k,n-k);
    }
    return t;
}

void join(node* t) // слияние узлов t и t->link
{
    node* p = t->link;
    int* a = new int[t->len+p->len];
    std::copy(a, a + t->len, t->key);
    std::copy(a+t->len,a+t->len + p->len, p->key);
    delete[] t->key;
    t->key = a;
    t->len += p->len;
    t->link = p->link;
    delete p;
}


template <typename T>
bool pdf(std::vector<T> x)
{
    T rez1 = 0, rez2 = 0;
    for(auto i : x)
    {
        rez1 += std::pow(i+1, 2.0);
        rez2 += std::pow(i-1, 2.0);
    }
    return (std::sqrt(rez1) < 1.5 || std::sqrt(rez2) < 1.5) && std::sqrt(rez2) > 0.5 ? true : false;
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
        struct TrieNode *visited,
        std::vector<std::vector<double>> &samples,
        std::vector<double> dx,
        size_t &counter,
        size_t &fe_count)
{
    while(!points.empty())
    {
        auto t = points.back();
        points.pop_back();

        if(search(visited,t))
        {
            counter++;
            continue;
        }
        insert(visited,t);

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


void FloodFill_MultipleGrids_VonNeumann_trie2(std::vector<std::vector<double>>& grids,
        std::vector<std::vector<int>> &points,
        node *visited,
        std::vector<std::vector<double>> &samples,
        std::vector<double> dx,
        size_t &counter,
        size_t &fe_count)
{
    int *tt = new int [dimension];

    while(!points.empty())
    {
        auto t = points.back();
        points.pop_back();

        for(size_t i = 0; i != dimension; i++)
        {
            tt[i] = t[i];
        }
        if(find(visited,tt) == 0)
        {
            counter++;
            continue;
        }
        visited = insert(visited,tt);

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

    delete [] tt;
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

    TrieNode *visited = getNode();

    std::vector<std::vector<double> > samples;

    size_t counter = 0;
    size_t fe_count = 0;

    FloodFill_MultipleGrids_VonNeumann_trie(grids, points, visited, samples, dx, counter, fe_count);

    std::cout << counter << std::endl;
    std::cout << "fe count: " << fe_count << std::endl;
    std::cout << "samples: " << samples.size() << std::endl;
    std::cout << samples.size()/double(fe_count) << std::endl;

    //print2file2d("maps/sample2d.dat", samples);
}


void b4MultipleGrids_trie2(std::vector<double> init_point)
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

    node *visited;

    std::vector<std::vector<double> > samples;

    size_t counter = 0;
    size_t fe_count = 0;

    FloodFill_MultipleGrids_VonNeumann_trie2(grids, points, visited, samples, dx, counter, fe_count);

    std::cout << counter << std::endl;
    std::cout << "fe count: " << fe_count << std::endl;
    std::cout << "samples: " << samples.size() << std::endl;
    std::cout << samples.size()/double(fe_count) << std::endl;

    //print2file2d("maps/sample2d.dat", samples);
}


int main()
{
    size_t dim = dimension;
    std::vector<double> start(dim, -1);

    Timer time_cpp11;

    time_cpp11.reset();
    b4MultipleGrids(start);
    std::cout << time_cpp11.elapsed_seconds() << std::endl;
    time_cpp11.reset();
    b4MultipleGrids_trie(start);
    std::cout << time_cpp11.elapsed_seconds() << std::endl;
    time_cpp11.reset();
    b4MultipleGrids_trie2(start);
    std::cout << time_cpp11.elapsed_seconds() << std::endl;
}
