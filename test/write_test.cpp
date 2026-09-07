#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/write/write.hpp"
#include <cstdio>
#include <doctest/doctest.h>
#include <highfive/highfive.hpp>
#include <string>

TEST_CASE("testing save_sim") {
  const std::string path = "/tmp/femib_write_test_sim.h5";

  femib::write::save_sim(path, "test_sim");

  HighFive::File file(path, HighFive::File::ReadOnly);
  std::string sim_name;
  file.getAttribute("sim_name").read(sim_name);
  CHECK(sim_name == "test_sim");

  std::remove(path.c_str());
}

TEST_CASE("testing save_plot_data") {
  const std::string path = "/tmp/femib_write_test_plot_data.h5";
  femib::write::save_sim(path, "test_sim");

  femib::write::plot_data<float, 2> data0;
  data0.time = 0;
  data0.x = {{0.0f, 0.0f}, {1.0f, 1.0f}};
  data0.u = {{1.0f, 2.0f}, {3.0f, 4.0f}};
  // q empty
  femib::write::save_plot_data(path, data0);

  femib::write::plot_data<float, 2> data1;
  data1.time = 1;
  data1.x = {{0.5f, 0.5f}};
  data1.u = {{5.0f, 6.0f}};
  data1.q = {{Eigen::Matrix<float, 1, 1>(7.0f)}};
  femib::write::save_plot_data(path, data1);

  HighFive::File file(path, HighFive::File::ReadOnly);

  std::vector<std::vector<float>> x_read;
  file.getDataSet("timestep_0/x").read(x_read);
  CHECK(x_read.size() == 2);
  CHECK(x_read[0][0] == doctest::Approx(0.0));
  CHECK(x_read[1][1] == doctest::Approx(1.0));

  std::vector<std::vector<float>> u_read;
  file.getDataSet("timestep_0/u").read(u_read);
  CHECK(u_read.size() == 2);
  CHECK(u_read[0][0] == doctest::Approx(1.0));
  CHECK(u_read[0][1] == doctest::Approx(2.0));
  CHECK(u_read[1][0] == doctest::Approx(3.0));
  CHECK(u_read[1][1] == doctest::Approx(4.0));

  CHECK_FALSE(file.exist("timestep_0/q"));

  std::vector<std::vector<float>> q_read;
  file.getDataSet("timestep_1/q").read(q_read);
  CHECK(q_read.size() == 1);
  CHECK(q_read[0][0] == doctest::Approx(7.0));

  std::vector<std::vector<float>> x1_read;
  file.getDataSet("timestep_1/x").read(x1_read);
  CHECK(x1_read.size() == 1);
  CHECK(x1_read[0][0] == doctest::Approx(0.5));

  std::string sim_name;
  file.getAttribute("sim_name").read(sim_name);
  CHECK(sim_name == "test_sim");

  std::remove(path.c_str());
}
