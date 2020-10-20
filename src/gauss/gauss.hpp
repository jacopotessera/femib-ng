#ifndef GAUSS_HPP_INCLUDED_
#define GAUSS_HPP_INCLUDED_

#include <vector>
#include <functional>

namespace femib::gauss {

	template<typename T=float, int d=2>
	struct dvec {
		T x;
		T y;
	};

	template<typename T, int d>
	struct node {
        	T weight;
        	dvec<T,d> node;
	};

	template<typename T, int d>
	struct rule {
        	std::vector<node<T,d>> nodes;
	};

	template<typename T, int d>
	T integrate(const rule<T,d>rule, std::function<T(dvec<T,d>)> &f);

}
#endif
