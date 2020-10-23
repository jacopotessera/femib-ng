#ifndef AFFINE_H_INCLUDED_
#define AFFINE_H_INCLUDED_

#include <functional>
#include <vector>

#include <Eigen/Core>

typedef Eigen::Matrix<float, 2, 2> dmat;
typedef Eigen::Matrix<float, 2, 1> dvec;
typedef std::vector<dvec> dtrian;

dvec affine(const dtrian &t, const dvec &x);
dvec affine_inv(const dtrian &t, const dvec &x);
dmat affineB(const dtrian &t);
dvec affineb(const dtrian &t);

#endif
