#include "gauss.hpp"

#include <functional>

template<typename T, int d>
T femib::gauss::integrate(const femib::gauss::rule<T,d> rule, std::function<T(femib::gauss::dvec<T,d>)> &f) { 
	T integral = 0;
        for(femib::gauss::node<T,d> node : rule.nodes) {
		integral += node.weight*f(node.node);
	}
	return integral;
}

template float femib::gauss::integrate<float,2>(const femib::gauss::rule<float,2> rule, std::function<float(femib::gauss::dvec<float,2>)> &f);
