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

document plot_data2doc(femib::mongo::plot_data t) {
  document data_builder{};
  data_builder << "id" << t.id << "time" << t.time;

  auto array_builder_u = data_builder << "u" << open_array;
  for (int i = 0; i < t.u.size(); ++i) {
    auto doc = array_builder_u << bsoncxx::builder::stream::open_array;
    for (int j = 0; j < t.u[i].size(); ++j) {
      auto doc_ = doc << bsoncxx::builder::stream::open_array;
      for (int k = 0; k < t.u[i][j].size(); ++k) {
        doc_ << t.u[i][j][k];
      }
      doc_ << bsoncxx::builder::stream::close_array;
    }
    doc << bsoncxx::builder::stream::close_array;
  }
  array_builder_u << close_array;

  auto array_builder_q = data_builder << "q" << open_array;
  for (int i = 0; i < t.q.size(); ++i) {
    auto doc = array_builder_q << bsoncxx::builder::stream::open_array;
    for (int j = 0; j < t.q[i].size(); ++j) {
      auto doc_ = doc << bsoncxx::builder::stream::open_array;
      for (int k = 0; k < t.q[i][j].size(); ++k) {
        doc_ << t.q[i][j][k];
      }
      doc_ << bsoncxx::builder::stream::close_array;
    }
    doc << bsoncxx::builder::stream::close_array;
  }
  array_builder_q << close_array;

  auto array_builder_x = data_builder << "x" << open_array;
  for (int i = 0; i < t.x.size(); ++i) {
    auto doc = array_builder_x << bsoncxx::builder::stream::open_array;
    for (int j = 0; j < t.x[i].size(); ++j) {
      auto doc_ = doc << bsoncxx::builder::stream::open_array;
      for (int k = 0; k < t.x[i][j].size(); ++k) {
        doc_ << t.x[i][j][k];
      }
      doc_ << bsoncxx::builder::stream::close_array;
    }
    doc << bsoncxx::builder::stream::close_array;
  }
  array_builder_x << close_array;

  return data_builder;
}

void femib::mongo::save_plot_data(std::string dbname,
                                  femib::mongo::plot_data t) {

  mongocxx::instance inst{};
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
