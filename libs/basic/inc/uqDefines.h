#ifndef __UQ_DEFINES_H__
#define __UQ_DEFINES_H__

#include <iostream>
#include <limits>

int uqMyRank();

#define UQ_UNAVAILABLE_RANK -1
#define UQ_INVALID_DOF_ID   UINT_MAX
#define UQ_INVALID_NODE_ID  UINT_MAX

#define UQ_OK_RC                          0
#define UQ_INCOMPLETE_IMPLEMENTATION_RC  -1
#define UQ_INVALID_INPUT_PARAMETER_RC    -2
#define UQ_INVALID_INTERNAL_RESULT_RC    -3
#define UQ_INVALID_INTERNAL_STATE_RC     -4
#define UQ_FAILED_TO_OPEN_FILE_RC        -5
#define UQ_MATRIX_IS_NOT_POS_DEFINITE_RC -6
#define UQ_FAILED_READING_FILE_RC        -7

#define UQ_RC_MACRO(iRC,givenRank,where,what,retValue) \
  if (iRC) {                                           \
    int macroRank = givenRank;                         \
    if (macroRank == UQ_UNAVAILABLE_RANK) {            \
      macroRank = uqMyRank();                          \
    }                                                  \
    std::cerr << "UQ RC ERROR"                         \
              << ", rank "  << macroRank               \
              << ", in "    << where                   \
              << ": "       << what                    \
              << ", iRC = " << iRC                     \
              << std::endl;                            \
    return retValue;                                   \
  }

#define UQ_TEST_MACRO(test,givenRank,where,what,retValue) \
  if (test) {                                             \
    int macroRank = givenRank;                            \
    if (macroRank == UQ_UNAVAILABLE_RANK) {               \
      macroRank = uqMyRank();                             \
    }                                                     \
    std::cerr << "UQ TEST ERROR"                          \
              << ", rank " << macroRank                   \
              << ", in "   << where                       \
              << ": "      << what                        \
              << std::endl;                               \
    return retValue;                                      \
  }

#define UQ_FATAL_RC_MACRO(iRC,givenRank,where,what) \
  if (iRC) {                                        \
    int macroRank = givenRank;                      \
    if (macroRank == UQ_UNAVAILABLE_RANK) {         \
      macroRank = uqMyRank();                       \
    }                                               \
    std::cerr << "UQ RC FATAL ERROR"                \
              << ", rank "  << macroRank            \
              << ", in "    << where                \
              << ": "       << what                 \
              << ", iRC = " << iRC                  \
              << ". Exiting..."                     \
              << std::endl;                         \
    exit(1);                                        \
  }

#define UQ_FATAL_TEST_MACRO(test,givenRank,where,what) \
  if (test) {                                          \
    int macroRank = givenRank;                         \
    if (macroRank == UQ_UNAVAILABLE_RANK) {            \
      macroRank = uqMyRank();                          \
    }                                                  \
    std::cerr << "UQ TEST FATAL ERROR"                 \
              << ", rank "  << macroRank               \
              << ", in "    << where                   \
              << ": "       << what                    \
              << ". Exiting..."                        \
              << std::endl;                            \
    exit(1);                                           \
  }

#endif // __UQ_DEFINES_H__
