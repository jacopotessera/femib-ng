/*
*	affine.cpp
*/

#include "affine.hpp"

dmat affineB(const dtrian &t) {
    dmat B;
    B << 	
        t[1](0) - t[0](0), t[2](0) - t[0](0), 
	t[1](1) - t[0](1), t[2](1) - t[0](1);
    return B;
}
dvec affineb(const dtrian &t) {
    dvec b = t[0];
    return b;
}

dvec affine(const dtrian &t, const dvec &x) {
    dvec y = affineB(t)*x + affineb(t);
    return y;
}

dvec affine_inv(const dtrian &t, const dvec &x) {
    dvec y = affineB(t).inverse()*(x-affineb(t));
    return y;
}

std::function<dvec(const dvec&)> affine(const dtrian &t) {
	return [&t](const dvec &x){
		return affine(t, x);
	};
}

/*std::function<dvec(const dvec&)> affine_inv(const dtrian &t) {
	
	return [&t](const dvec &x){
		dmat B;
		B << 	t[1](0) - t[0](0), t[2](0) - t[0](0),
		   	t[1](1) - t[0](1), t[2](1) - t[0](1);
		dvec b = t[0];
		return B*x+b;
	};
}*/
