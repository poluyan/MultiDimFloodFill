#ifndef MOORE_H
#define MOORE_H


template <typename T>
bool increase(const std::vector<std::vector<T>>& v, std::vector<size_t>& it)
{
    for(size_t i = 0, size = it.size(); i != size; i++)
    {
        const size_t index = size - 1 - i;
        ++it[index];
        if(it[index] == v[index].size())
        {
            it[index] = 0;
        }
        else
        {
            return true;
        }
    }
    return false;
}

template <typename T>
std::vector<T> get_line(const std::vector<std::vector<T>>& v, std::vector<size_t>& it)
{
    std::vector<T> rez(v.size());
    for(size_t i = 0, size = v.size(); i != size; i++)
    {
        rez[i] = v[i][it[i]];
    }
    return rez;
}

template <typename T>
std::vector<std::vector<T>> iterate(const std::vector<std::vector<T>>& v)
{
    std::vector<size_t> it(v.size(), 0);
    std::vector<std::vector<T>> values;
    do
    {
        values.push_back(do_job(v, it));
    }
    while(increase(v, it));
    return values;
}

#endif