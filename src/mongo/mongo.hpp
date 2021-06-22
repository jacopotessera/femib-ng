#ifndef MONGO_HPP_INCLUDED_
#define MONGO_HPP_INCLUDED_

#include <string>
#include <vector>

namespace femib::mongo {

struct plot_data {
  std::string id;
  int time;
  std::vector<std::vector<std::vector<double>>> u;
  std::vector<std::vector<std::vector<double>>> q;
  std::vector<std::vector<std::vector<double>>> x;
};

void save_plot_data(std::string dbname, plot_data t);

} // namespace femib::mongo
#endif
