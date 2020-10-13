#include "affine/affine.hpp"
#include <iostream>

int main(){
    std::cout << "Hello CMake!" << std::endl;

    dvec p1;
    p1 << 1,2;
    dvec p2;
    p2 << 3,0;
    dvec p3;
    p3 << 2,4;

    dvec x;
    x << 1,0;
    
    std::vector<dvec> t;
    t.push_back(p1);
    t.push_back(p2);
    t.push_back(p3);

    dvec y = affine(t, x);
    std::cout << y << std::endl;

    dvec z = affine_inv(t, y);
    std::cout << z << std::endl;

    return 0;
}
