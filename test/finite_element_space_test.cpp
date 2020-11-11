#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <Eigen/Sparse>
#include <spdlog/spdlog.h>
#include "../src/finite_element_space/finite_element_space.hpp"
#include "../src/gauss/gauss.hpp"
#include "../src/mesh/mesh.hpp"
#include "../src/affine/affine.hpp"

int get_index(int i, int n)
{
  if(n == 0){
    if(i == 0){
      return 0;
    } else if(i==1){
      return 1;
    } else if(i==2){
      return 4;
    }
  } else if(n == 1){
     if(i == 0){
       return 1;
     } else if(i==1){
       return 2;
     } else if(i==2){
        return 4;
     }
   } else if(n == 2){
     if(i == 0){
       return 2;
     } else if(i==1){
       return 3;
     } else if(i==2){
       return 4;
     }
   } else if(n == 3){
      if(i == 0){
        return 3;
      } else if(i==1){
        return 0;
      } else if(i==2){
        return 4;
      }
    }
  return -1;
}

TEST_CASE("testing finite_element_space") {

  std::vector<Eigen::Triplet<float>> M;

  femib::gauss::node<float, 2> node1 = {1.0 / 6.0, {1.0 / 6.0, 1.0 / 6.0}};
  femib::gauss::node<float, 2> node2 = {1.0 / 6.0, {1.0 / 6.0, 2.0 / 3.0}};
  femib::gauss::node<float, 2> node3 = {1.0 / 6.0, {2.0 / 3.0, 1.0 / 6.0}};
  femib::gauss::rule<float, 2> rule = {{node1, node2, node3}};

  femib::types::mesh<float, 2> mesh = {
      // P
      {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.5, 0.5}},
      // T
      {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}

      },
      // E
      {0, 1, 2, 3}};

  femib::finite_element::finite_element<float, 2, 1> f =
      femib::finite_element::create_finite_element_P1_2d1d<float, 2, 1>();

  femib::finite_element_space::finite_element_space<float, 2, 1> s = {f, mesh};

	for(int n=0;n<s.mesh.T.size();++n)
	{
		for(int i=0;i<s.finite_element.base_functions.size();++i)
		{
			for(int j=0;j<s.finite_element.base_functions.size();++j)
			{
					femib::types::F<float,2,1> a;
					femib::types::F<float,2,1> b;
					femib::types::dtrian<float,2> t = s.mesh[n];
					a.dx = [&](const femib::types::dvec<float,2> &x){
						return (affineBinv(t)*f.base_functions[i].dx(affineBinv(t)*(x-affineb(t))));};
					b.dx = [&](const femib::types::dvec<float,2> &x){
						return (affineBinv(t)*f.base_functions[j].dx(affineBinv(t)*(x-affineb(t))));};
				//float m = femib::mesh::integrate<float,2>(rule, [](femib::types::dvec<float,2>){return (float)1.0;}, t);
				float m = femib::mesh::integrate<float,2>(rule, [&](femib::types::dvec<float,2> x){
						float a_0 = (affineBinv(t)*f.base_functions[i].dx(affineBinv(t)*(x-affineb(t))))[0];
						float a_1 = (affineBinv(t)*f.base_functions[i].dx(affineBinv(t)*(x-affineb(t))))[1];
						float b_0 = (affineBinv(t)*f.base_functions[j].dx(affineBinv(t)*(x-affineb(t))))[0];
						float b_1 = (affineBinv(t)*f.base_functions[j].dx(affineBinv(t)*(x-affineb(t))))[1];
						return a_0*b_0 + a_1*b_1;
						//return a.dx(x)[0]*b.dx(x)[0] + a.dx(x)[1]*b.dx(x)[1];
						}, t);
				std::cout << "n: " << n << ", i: " << i << ", j: " << j << " -> " << m << std::endl;
				M.push_back(Eigen::Triplet<float>(get_index(i,n),get_index(j,n),m));
			}
		}
	}

	Eigen::SparseMatrix<float> sM = Eigen::SparseMatrix<float>(5,5);
	sM.setFromTriplets(M.begin(),M.end());
	std::cout << sM << std::endl;

    CHECK(true);
}
