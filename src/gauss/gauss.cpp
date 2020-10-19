#include <vector>
#include <functional>

template<typename T=float, int d=2>
struct dvec {
	
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
T integrate(const rule<T,d> rule, std::function<T(dvec<T,d>)> &f) { 
	T integral = 0;
        for(node<T,d> node : rule) {
		integral += node.weight*f(node.node);
	}
	return integral;
}

