#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/mongo/mongo.hpp"
#include <doctest/doctest.h>
#include <string>

TEST_CASE("testing mongo") {

  std::vector<std::vector<std::vector<float>>> u = {
      {{0, 0}, {0, 0}},   {{0.5, 0}, {1, 1}},   {{1, 0}, {0, 0}},
      {{0, 0.5}, {1, 1}}, {{0.5, 0.5}, {2, 2}}, {{1, 0.5}, {1, 1}},
      {{0, 1}, {0, 0}},   {{0.5, 1}, {1, 1}},   {{1, 1}, {0, 0}}};
  std::vector<std::vector<std::vector<float>>> q;
  std::vector<std::vector<std::vector<float>>> x;

  std::string dbname = "femib_test";
  std::string id = "666";
  int time = 0;
  femib::mongo::save_sim(dbname, id);
  femib::mongo::plot_data p = {id, time, u, q, x};
  femib::mongo::save_plot_data(dbname, p);
  CHECK(true);
}
