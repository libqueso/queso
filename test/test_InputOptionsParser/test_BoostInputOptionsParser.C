#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
#include <queso/asserts.h>
#include <queso/BoostInputOptionsParser.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <set>
#include <vector>

#define TEST_INT -1
#define TEST_UINT 1
#define TEST_DOUBLE -1.0
#define TEST_BOOL true
#define TEST_STRING "hello world"
#define TEST_SET_UINT "1 2 3"
#define TEST_VECTOR_DOUBLE "-1.0 2.0 3.0"

class TestOptionsValues
{
public:
  TestOptionsValues(const std::string & filename);

  QUESO::BoostInputOptionsParser * m_parser;

  std::string m_option_int;
  std::string m_option_unsigned_int;
  std::string m_option_double;
  std::string m_option_bool;
  std::string m_option_string;
  std::string m_option_set_unsigned_int;
  std::string m_option_vector_double;

  int m_int;
  unsigned int m_unsigned_int;
  double m_double;
  bool m_bool;
  std::string m_string;
  std::set<unsigned int> m_set_unsigned_int;
  std::vector<double> m_vector_double;
};

TestOptionsValues::TestOptionsValues(const std::string & filename)
  :
    m_parser(new QUESO::BoostInputOptionsParser(filename)),
    m_option_int("int"),
    m_option_unsigned_int("unsigned_int"),
    m_option_double("double"),
    m_option_bool("bool"),
    m_option_string("std_string"),
    m_option_set_unsigned_int("set_unsigned_int"),
    m_option_vector_double("vector_double")
{
  // Register all options with parser
  m_parser->registerOption<int>(m_option_int, TEST_INT, "integer");
  m_parser->registerOption<unsigned int>(m_option_unsigned_int, TEST_UINT,
      "unsigned integer");
  m_parser->registerOption<double>(m_option_double, TEST_DOUBLE, "double");
  m_parser->registerOption<bool>(m_option_bool, TEST_BOOL, "boolean");
  m_parser->registerOption<std::string>(m_option_string, TEST_STRING,
      "string");
  m_parser->registerOption<std::string>(m_option_set_unsigned_int,
      TEST_SET_UINT, "set of unsigned integers");
  m_parser->registerOption<std::string>(m_option_vector_double,
      TEST_VECTOR_DOUBLE, "vector of doubles");

  // Read the input file
  m_parser->scanInputFile();

  m_parser->getOption<int>(m_option_int, m_int);
  m_parser->getOption<unsigned int>(m_option_unsigned_int, m_unsigned_int);
  m_parser->getOption<double>(m_option_double, m_double);
  m_parser->getOption<bool>(m_option_bool, m_bool);
  m_parser->getOption<std::string>(m_option_string, m_string);
  m_parser->getOption<std::set<unsigned int> >(m_option_set_unsigned_int,
      m_set_unsigned_int);
  m_parser->getOption<std::vector<double> >(m_option_vector_double,
      m_vector_double);
}

void test_default()
{
  std::string inputFileName =
    "test_InputOptionsParser/test_options_default.txt";

  // Find correct file for out-of-source builds
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    inputFileName = test_srcdir + ('/' + inputFileName);

  TestOptionsValues optionsValues(inputFileName);

  queso_require_equal_to(optionsValues.m_int, -1);
  queso_require_equal_to(optionsValues.m_unsigned_int, 1);
  queso_require_equal_to(optionsValues.m_double, -1.0);
  queso_require_equal_to(optionsValues.m_bool, true);
  queso_require_equal_to(optionsValues.m_string, std::string("hello world"));
  queso_require_equal_to(optionsValues.m_set_unsigned_int.size(), 3);
  queso_require_equal_to(optionsValues.m_set_unsigned_int.count(1), 1);
  queso_require_equal_to(optionsValues.m_set_unsigned_int.count(2), 1);
  queso_require_equal_to(optionsValues.m_set_unsigned_int.count(3), 1);
  queso_require_equal_to(optionsValues.m_vector_double.size(), 3);
  queso_require_equal_to(optionsValues.m_vector_double[0], -1.0);
  queso_require_equal_to(optionsValues.m_vector_double[1], 2.0);
  queso_require_equal_to(optionsValues.m_vector_double[2], 3.0);
}

void test_good()
{
  std::string inputFileName =
    "test_InputOptionsParser/test_options_good.txt";

  // Find correct file for out-of-source builds
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    inputFileName = test_srcdir + ('/' + inputFileName);

  TestOptionsValues optionsValues(inputFileName);

  queso_require_equal_to(optionsValues.m_int, -2);
  queso_require_equal_to(optionsValues.m_unsigned_int, 2);
  queso_require_equal_to(optionsValues.m_double, -2.0);
  queso_require_equal_to(optionsValues.m_bool, false);
  queso_require_equal_to(optionsValues.m_string, std::string("hello"));
  queso_require_equal_to(optionsValues.m_set_unsigned_int.size(), 3);
  queso_require_equal_to(optionsValues.m_set_unsigned_int.count(2), 1);
  queso_require_equal_to(optionsValues.m_set_unsigned_int.count(3), 1);
  queso_require_equal_to(optionsValues.m_set_unsigned_int.count(4), 1);
  queso_require_equal_to(optionsValues.m_vector_double.size(), 3);
  queso_require_equal_to(optionsValues.m_vector_double[0], -2.0);
  queso_require_equal_to(optionsValues.m_vector_double[1], 3.0);
  queso_require_equal_to(optionsValues.m_vector_double[2], 4.0);
}

void test_empty()
{
  std::string inputFileName =
    "test_InputOptionsParser/test_options_bad.txt";

  // Find correct file for out-of-source builds
  const char * test_srcdir = std::getenv("srcdir");
  if (test_srcdir)
    inputFileName = test_srcdir + ('/' + inputFileName);

  TestOptionsValues optionsValues(inputFileName);

  queso_require_equal_to(optionsValues.m_set_unsigned_int.empty(), true);
  queso_require_equal_to(optionsValues.m_vector_double.empty(), true);
}
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS

int main(int argc, char ** argv)
{
#ifndef DISABLE_BOOST_PROGRAM_OPTIONS
  test_default();
  test_good();
  test_empty();

  return 0;
#else
  return 77;
#endif  // DISABLE_BOOST_PROGRAM_OPTIONS
}
