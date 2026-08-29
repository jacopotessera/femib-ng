#ifndef FEMIB_NAME_REPORTER_H
#define FEMIB_NAME_REPORTER_H

#include <doctest/doctest.h>
#include <iostream>

using namespace doctest;

struct NameReporter : IReporter {
  std::ostream &stdout_stream;
  const ContextOptions &opt;
  int index = 1;

  NameReporter(const ContextOptions &in) : stdout_stream(*in.cout), opt(in) {}

  void test_case_start(const TestCaseData &in) override {
    // TODO use fmt/color.h
    stdout_stream << "\033[1;36m" << "[doctest] " << "\033[0m" << index << ": "
                  << in.m_name << '\n';
    index++;
  }

  void test_case_end(const CurrentTestCaseStats &) override {}
  void report_query(const QueryData &) override {}
  void test_run_start() override {}
  void test_run_end(const TestRunStats &) override {}
  void test_case_reenter(const TestCaseData &) override {}
  void log_assert(const AssertData &) override {}
  void log_message(const MessageData &) override {}
  void test_case_exception(const TestCaseException &) override {}
  void subcase_start(const SubcaseSignature &) override {}
  void subcase_end() override {}
  void test_case_skipped(const TestCaseData &) override {}
};

REGISTER_REPORTER("names", 1, NameReporter);

#endif // FEMIB_NAME_REPORTER_H
