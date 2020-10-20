#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "affine.hpp"

TEST_CASE("testing affine") {
	dvec vertex_1 = {1.0, 2.0};
	dvec vertex_2 = {3.0, 0.0};
	dvec vertex_3 = {2.0, 4.0};
	std::vector<dvec> triangle = {vertex_1, vertex_2, vertex_3};

	dvec point_x = {1.0, 0.0};
	dvec point_y = affine(triangle, point_x);
	CHECK(point_y(0) == 3);
	CHECK(point_y(1) == 0);

	dvec point_z = affine_inv(triangle, point_y);
	CHECK(point_x == point_z);
}

