#ifndef PDF_H
#define PDF_H

template <typename T>
bool pdf(std::vector<T> x)
{
//    T rez1 = 0, rez2 = 0;
//    for(auto i : x)
//    {
//        rez1 += std::pow(i+0.5, 2.0);
//        rez2 += std::pow(i-0.5, 2.0);
//    }
//    return (std::sqrt(rez1) < 1.5 || std::sqrt(rez2) < 1.5) && std::sqrt(rez2) > 0.5 ? true : false;

//    T sum1 = 0, sum2 = 0;
//    for(auto i : x)
//    {
//        sum1 += i*i;
//        sum2 += (i+1)*(i+1);
//    }
//
//    return std::exp(-5.0*std::abs(sum1 - 1.0)) > 0.1 && std::sqrt(sum2) > 1.75 ? true : false ;

    T rez = 0;
    if(x.size() == 2)
    {
        T a = 2.0, b = 0.5;
        rez = std::signbit(std::pow(x[0]*x[0]+x[1]*x[1]-(a*a+b*b),2.0)-4.0*a*a*(b*b));
    }
    if(x.size() == 3)
    {
        T a = 2.0, b = 0.5;
        rez = std::signbit(std::pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]-(a*a+b*b),2.0)-4.0*a*a*(b*b-x[2]*x[2]));
    }
    return rez;
}

#endif