#include <benchmark/benchmark.h>

#include "../src/gauss/gauss.hpp"
#include <assert.h>
#include <cmath>

const float EPSILON = std::numeric_limits<float>::epsilon();

const femib::gauss::node<float, 2> node1 = {1.0 / 6.0, {1.0 / 6.0, 1.0 / 6.0}};
const femib::gauss::node<float, 2> node2 = {1.0 / 6.0, {1.0 / 6.0, 2.0 / 3.0}};
const femib::gauss::node<float, 2> node3 = {1.0 / 6.0, {2.0 / 3.0, 1.0 / 6.0}};
const femib::gauss::rule<float, 2> rule = {{node1, node2, node3}};
std::function<float(femib::types::dvec<float, 2>)> f =
    [](femib::types::dvec<float, 2> x) { return x(0) + x(1); };

static void gauss_integrate(benchmark::State &state) {
  for (auto _ : state) {
    float integral = femib::gauss::integrate<float, 2>(rule, f);
    assert(std::fabs(integral - 1.0 / 3.0) <= EPSILON);
  }
}

BENCHMARK(gauss_integrate);

BENCHMARK_MAIN();
