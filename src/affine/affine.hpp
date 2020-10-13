/*
*	affine.h
*/

#ifndef AFFINE_H_INCLUDED_
#define AFFINE_H_INCLUDED_

#include <functional>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

typedef Eigen::Matrix<float, 2, 2> dmat;
typedef Eigen::Matrix<float, 2, 1> dvec;
typedef std::vector<dvec> dtrian;

dvec affine(const dtrian &t, const dvec &x);
dmat affineB(const dtrian &t);
dvec affineb(const dtrian &t);


dvec affine_inv(const dtrian &t, const dvec &x);

std::function<dvec(const dvec&)> affine(const dtrian &t);
//std::function<dvec(const dvec&)> affine_inv(const dtrian &t);

#endif

