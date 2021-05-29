#include "../src/affine/affine.hpp"
#include "../src/femib/femib.hpp"
#include "../src/finite_element/P1_2d1d.hpp"
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/gauss/gauss_lagrange_2_2d.hpp"
#include "../src/mesh/mesh.hpp"
#include "../src/read/read.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace femib::affine;

int get_index(const femib::types::nodes<float, 2> &nodes, int i, int n) {
  return nodes.T[n][i];
}

int main() {

  femib::gauss::rule<float, 2> rule =
      femib::gauss::create_gauss_2_2d<float, 2>();
  std::string mesh_dir = MESH_DIR;
  femib::types::mesh<float, 2> mesh = femib::mesh::read<float, 2>(
      mesh_dir + "p5.mat", mesh_dir + "t5.mat", mesh_dir + "e5.mat");

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s = {f, mesh};
  s.nodes = f.build_nodes(mesh);

  femib::poisson::poisson<float, 2, 1> poisson;
  poisson.V = s;

  poisson.M =
      [&f, &rule](
          femib::finite_element_space::finite_element_space<float, 2, 1> s,
          femib::finite_element_space::finite_element_space<float, 2, 1> ss) {
        std::vector<Eigen::Triplet<float>> MM;
        for (int n = 0; n < s.mesh.T.size(); ++n) {
          for (int i = 0; i < s.finite_element.base_functions.size(); ++i) {
            for (int j = 0; j < s.finite_element.base_functions.size(); ++j) {
              femib::types::F<float, 2, 1> a;
              femib::types::F<float, 2, 1> b;
              femib::types::dtrian<float, 2> t = s.mesh[n];
              a.dx = [&](const femib::types::dvec<float, 2> &x) {
                return (affineBinv(t) * f.base_functions[i].dx(
                                            affineBinv(t) * (x - affineb(t))));
              };
              b.dx = [&](const femib::types::dvec<float, 2> &x) {
                return (affineBinv(t) * f.base_functions[j].dx(
                                            affineBinv(t) * (x - affineb(t))));
              };
              // float m = femib::mesh::integrate<float,2>(rule,
              // [](femib::types::dvec<float,2>){return (float)1.0;}, t);
              float m = femib::mesh::integrate<float, 2>(
                  rule,
                  [&t, &f, i, j](femib::types::dvec<float, 2> x) {
                    float a_0 = (affineBinv(t) *
                                 f.base_functions[i].dx(affineBinv(t) *
                                                        (x - affineb(t))))[0];
                    float a_1 = (affineBinv(t) *
                                 f.base_functions[i].dx(affineBinv(t) *
                                                        (x - affineb(t))))[1];
                    float b_0 = (affineBinv(t) *
                                 f.base_functions[j].dx(affineBinv(t) *
                                                        (x - affineb(t))))[0];
                    float b_1 = (affineBinv(t) *
                                 f.base_functions[j].dx(affineBinv(t) *
                                                        (x - affineb(t))))[1];
                    return a_0 * b_0 + a_1 * b_1;
                    // return a.dx(x)[0]*b.dx(x)[0] + a.dx(x)[1]*b.dx(x)[1];
                  },
                  t);
              float f_ = femib::mesh::integrate<float, 2>(
                  rule,
                  [&t, &f, i](femib::types::dvec<float, 2> x) {
                    float a_0 = -40 * x(0) * x(1) *
                                (f.base_functions[i].x(affineBinv(t) *
                                                       (x - affineb(t))))[0];
                    return a_0;
                  },
                  t);
              // std::cout << "n: " << n << ", i: " << i << ", j: " << j
              //          << " -> m: " << m << std::endl;
              // std::cout << "n: " << n << ", i: " << i << ", j: " << j
              //          << " -> f: " << f_ << std::endl;
              MM.push_back(Eigen::Triplet<float>(get_index(s.nodes, i, n),
                                                 get_index(s.nodes, j, n), m));
              // F.push_back(Eigen::Triplet<float>(get_index(s.nodes, i, n), 0,
              // f_));
            }
          }
        }
        return MM;
      };

  poisson.f = [&f, &rule](
                  femib::finite_element_space::finite_element_space<float, 2, 1>
                      s) {
    std::vector<Eigen::Triplet<float>> F;
    for (int n = 0; n < s.mesh.T.size(); ++n) {
      for (int i = 0; i < s.finite_element.base_functions.size(); ++i) {
        for (int j = 0; j < s.finite_element.base_functions.size(); ++j) {
          femib::types::F<float, 2, 1> a;
          femib::types::F<float, 2, 1> b;
          femib::types::dtrian<float, 2> t = s.mesh[n];
          a.dx = [&](const femib::types::dvec<float, 2> &x) {
            return (affineBinv(t) *
                    f.base_functions[i].dx(affineBinv(t) * (x - affineb(t))));
          };
          b.dx = [&](const femib::types::dvec<float, 2> &x) {
            return (affineBinv(t) *
                    f.base_functions[j].dx(affineBinv(t) * (x - affineb(t))));
          };
          // float m = femib::mesh::integrate<float,2>(rule,
          // [](femib::types::dvec<float,2>){return (float)1.0;}, t);
          float m = femib::mesh::integrate<float, 2>(
              rule,
              [&t, &f, i, j](femib::types::dvec<float, 2> x) {
                float a_0 =
                    (affineBinv(t) * f.base_functions[i].dx(
                                         affineBinv(t) * (x - affineb(t))))[0];
                float a_1 =
                    (affineBinv(t) * f.base_functions[i].dx(
                                         affineBinv(t) * (x - affineb(t))))[1];
                float b_0 =
                    (affineBinv(t) * f.base_functions[j].dx(
                                         affineBinv(t) * (x - affineb(t))))[0];
                float b_1 =
                    (affineBinv(t) * f.base_functions[j].dx(
                                         affineBinv(t) * (x - affineb(t))))[1];
                return a_0 * b_0 + a_1 * b_1;
                // return a.dx(x)[0]*b.dx(x)[0] + a.dx(x)[1]*b.dx(x)[1];
              },
              t);
          float f_ = femib::mesh::integrate<float, 2>(
              rule,
              [&t, &f, i](femib::types::dvec<float, 2> x) {
                float a_0 = -40 * x(0) * x(1) *
                            (f.base_functions[i].x(affineBinv(t) *
                                                   (x - affineb(t))))[0];
                return a_0;
              },
              t);
          // std::cout << "n: " << n << ", i: " << i << ", j: " << j
          //          << " -> m: " << m << std::endl;
          // std::cout << "n: " << n << ", i: " << i << ", j: " << j
          //          << " -> f: " << f_ << std::endl;
          // MM.push_back(Eigen::Triplet<float>(get_index(s.nodes, i, n),
          //                                  get_index(s.nodes, j, n), m));
          F.push_back(Eigen::Triplet<float>(get_index(s.nodes, i, n), 0, f_));
        }
      }
    }
    return F;
  };
}
