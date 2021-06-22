#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/mongo/mongo.hpp"
#include <doctest/doctest.h>

TEST_CASE("testing mongo") {

  std::vector<std::vector<std::vector<double>>> u = {
      {{0, 0}, {0, 0}},   {{0.5, 0}, {1, 1}},   {{1, 0}, {0, 0}},
      {{0, 0.5}, {1, 1}}, {{0.5, 0.5}, {2, 2}}, {{1, 0.5}, {1, 1}},
      {{0, 1}, {0, 0}},   {{0.5, 1}, {1, 1}},   {{1, 1}, {0, 0}}};
  std::vector<std::vector<std::vector<double>>> q;
  std::vector<std::vector<std::vector<double>>> x;

  femib::mongo::plot_data p = {"12345", 4, u, q, x};
  femib::mongo::save_plot_data("femib_test", p);
  CHECK(true);
}
