#include <iostream>
#include <algorithm>
#include <vector>

template <typename T>
T pdf(std::vector<T> x)
{
	T rez = 0;
	for(auto i : x)
	{
		rez += std::pow(i, 2.0);
	}
	return rez;
}

int main()
{
	
}