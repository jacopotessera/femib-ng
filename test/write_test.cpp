#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../src/write/write.hpp"
#include <cstdio>
#include <doctest/doctest.h>
#include <highfive/highfive.hpp>
#include <string>

TEST_CASE("testing save_sim") {
  const std::string path = "/tmp/femib_write_test_sim.h5";
  std::remove(path.c_str());
  femib::write::save_sim(path, "test_sim");

  HighFive::File file(path, HighFive::File::ReadOnly);
  std::string sim_name;
  file.getAttribute("sim_name").read(sim_name);
  CHECK(sim_name == "test_sim");
  std::string sim_type;
  file.getAttribute("sim_type").read(sim_type);
  CHECK(sim_type == "");

  std::remove(path.c_str());
}

TEST_CASE("testing save_sim with an explicit sim_type") {
  const std::string path = "/tmp/femib_write_test_sim_sim_type.h5";

  femib::write::save_sim(path, "test_sim", "ring",
                         femib::write::mode::overwrite);

  HighFive::File file(path, HighFive::File::ReadOnly);
  std::string sim_type;
  file.getAttribute("sim_type").read(sim_type);
  CHECK(sim_type == "ring");

  std::remove(path.c_str());
}

TEST_CASE("save_sim mode::create fails if the file already exists") {
  const std::string path = "/tmp/femib_write_test_create_guard.h5";
  std::remove(path.c_str());

  femib::write::save_sim(path, "test_sim");

  CHECK_THROWS(femib::write::save_sim(path, "test_sim"));

  std::remove(path.c_str());
}

TEST_CASE("save_sim mode::overwrite succeeds even if the file already exists") {
  const std::string path = "/tmp/femib_write_test_overwrite.h5";
  femib::write::save_sim(path, "first", "", femib::write::mode::overwrite);
  CHECK_NOTHROW(femib::write::save_sim(path, "second", "",
                                       femib::write::mode::overwrite));

  HighFive::File file(path, HighFive::File::ReadOnly);
  std::string sim_name;
  file.getAttribute("sim_name").read(sim_name);
  CHECK(sim_name == "second"); // proves the file was actually replaced

  std::remove(path.c_str());
}

TEST_CASE("save_sim mode::resume opens without truncating existing data") {
  const std::string path = "/tmp/femib_write_test_resume.h5";
  femib::write::save_sim(path, "test_sim", "", femib::write::mode::overwrite);

  femib::write::plot_data<float, 2> data0;
  data0.time = 0;
  data0.X = {{0.1f, 0.2f}};
  femib::write::save_plot_data(path, data0);

  // Resuming must NOT wipe out timestep_0, unlike mode::overwrite would.
  CHECK_NOTHROW(
      femib::write::save_sim(path, "test_sim", "", femib::write::mode::resume));

  HighFive::File file(path, HighFive::File::ReadOnly);
  CHECK(file.exist("timestep_0"));

  std::remove(path.c_str());
}

TEST_CASE("save_sim mode::resume fails if the file does not exist") {
  const std::string path = "/tmp/femib_write_test_resume_missing.h5";
  std::remove(path.c_str());
  CHECK_THROWS(
      femib::write::save_sim(path, "test_sim", "", femib::write::mode::resume));
}

TEST_CASE("testing save_plot_data") {
  const std::string path = "/tmp/femib_write_test_plot_data.h5";
  femib::write::save_sim(path, "test_sim", "", femib::write::mode::overwrite);

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

TEST_CASE("testing save_metadata") {
  const std::string path = "/tmp/femib_write_test_metadata.h5";
  femib::write::save_sim(path, "test_sim", "", femib::write::mode::overwrite);

  femib::write::save_metadata<double>(path, "deltat", 0.0002);
  femib::write::save_metadata<std::string>(path, "finite_element", "P1_B_2d2d");
  // overwrites if called again with the same key
  femib::write::save_metadata<double>(path, "deltat", 0.0005);

  HighFive::File file(path, HighFive::File::ReadOnly);
  HighFive::Group metadata = file.getGroup("metadata");

  double deltat_read;
  metadata.getAttribute("deltat").read(deltat_read);
  CHECK(deltat_read == doctest::Approx(0.0005));

  std::string fe_read;
  metadata.getAttribute("finite_element").read(fe_read);
  CHECK(fe_read == "P1_B_2d2d");

  std::remove(path.c_str());
}

TEST_CASE("testing save_mesh") {
  const std::string path = "/tmp/femib_write_test_mesh.h5";

  femib::types::mesh<double, 2> mesh;
  mesh.P = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}};
  mesh.T = {femib::types::ditrian<2>(0, 1, 2),
            femib::types::ditrian<2>(1, 3, 2)};
  mesh.E = {0, 1, 2, 3};

  femib::write::save_sim(path, "test_sim", "", femib::write::mode::overwrite);
  femib::write::save_mesh<double, 2>(path, mesh);

  HighFive::File file(path, HighFive::File::ReadOnly);

  std::vector<std::vector<double>> points_read;
  file.getDataSet("mesh/points").read(points_read);
  CHECK(points_read.size() == 4);
  CHECK(points_read[3][0] == doctest::Approx(1.0));
  CHECK(points_read[3][1] == doctest::Approx(1.0));

  std::vector<std::vector<int>> triangles_read;
  file.getDataSet("mesh/triangles").read(triangles_read);
  CHECK(triangles_read.size() == 2);
  CHECK(triangles_read[1][0] == 1);
  CHECK(triangles_read[1][1] == 3);
  CHECK(triangles_read[1][2] == 2);

  std::vector<int> edges_read;
  file.getDataSet("mesh/edges").read(edges_read);
  CHECK(edges_read.size() == 4);

  std::remove(path.c_str());
}
