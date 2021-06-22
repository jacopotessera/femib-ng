#include "spdlog/spdlog.h"

#include "mongo.hpp"

#include <bsoncxx/builder/stream/document.hpp>
#include <bsoncxx/json.hpp>
#include <bsoncxx/types.hpp>
#include <bsoncxx/types/value.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>

#include <bsoncxx/json.hpp>
#include <mongocxx/client.hpp>
#include <mongocxx/stdx.hpp>
#include <mongocxx/uri.hpp>

using bsoncxx::builder::stream::close_array;
using bsoncxx::builder::stream::close_document;
using bsoncxx::builder::stream::document;
using bsoncxx::builder::stream::finalize;
using bsoncxx::builder::stream::open_array;
using bsoncxx::builder::stream::open_document;

void add_array(document &data_builder, const std::string array_name,
               const std::vector<std::vector<std::vector<double>>> &u) {

  auto array_builder_u = data_builder << array_name << open_array;
  for (int i = 0; i < u.size(); ++i) {
    auto doc = array_builder_u << bsoncxx::builder::stream::open_array;
    for (int j = 0; j < u[i].size(); ++j) {
      auto doc_ = doc << bsoncxx::builder::stream::open_array;
      for (int k = 0; k < u[i][j].size(); ++k) {
        doc_ << u[i][j][k];
      }
      doc_ << bsoncxx::builder::stream::close_array;
    }
    doc << bsoncxx::builder::stream::close_array;
  }
  array_builder_u << close_array;
}

document plot_data2doc(const femib::mongo::plot_data &t) {
  document data_builder{};
  data_builder << "id" << t.id << "time" << t.time;

  add_array(data_builder, "u", t.u);
  add_array(data_builder, "q", t.q);
  add_array(data_builder, "x", t.x);

  return data_builder;
}

void femib::mongo::save_plot_data(std::string dbname,
                                  femib::mongo::plot_data t) {

  // mongocxx::instance inst{};
  mongocxx::client conn{mongocxx::uri{}};
  auto plot_data_collection = conn[dbname]["plot_data"];
  if (plot_data_collection.count_documents(
          document{} << "id" << t.id << "time" << t.time << finalize) == 0) {
    plot_data_collection.insert_one(plot_data2doc(t).view());
    SPDLOG_INFO("Saved plot_data id = {}, time = {}", t.id, t.time);
  } else {
    SPDLOG_ERROR("plot_data id = {}, time = {} already present in collection!",
                 t.id, t.time);
  }
}

void femib::mongo::save_sim(std::string dbname, std::string sim_name) {

  // mongocxx::instance inst{};
  mongocxx::client conn{mongocxx::uri{}};
  auto collection = conn[dbname]["sim"];
  if (collection.count_documents(document{} << "id" << sim_name << finalize) ==
      0) {
    collection.insert_one(document{} << "id" << sim_name << finalize);
    SPDLOG_INFO("Saved sim id = {}", sim_name);
  } else {
    SPDLOG_ERROR("sim id = {} already present in collection!", sim_name);
  }
}
