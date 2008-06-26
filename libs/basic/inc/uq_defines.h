#ifndef __UQ_DEFINES_H__
#define __UQ_DEFINES_H__

#include <iostream>
#include <limits>

#define UQ_UNAVAILABLE_RANK -1
#define UQ_INVALID_DOF_ID   UINT_MAX
#define UQ_INVALID_NODE_ID  UINT_MAX

#define UQ_OK_RC                         0
#define UQ_INCOMPLETE_IMPLEMENTATION_RC -1
#define UQ_INVALID_INPUT_PARAMETER_RC   -2
#define UQ_INVALID_INTERNAL_RESULT_RC   -3
#define UQ_INVALID_INTERNAL_STATE_RC    -4
#define UQ_FAILED_TO_OPEN_FILE_RC       -5
#define UQ_MATRIX_IS_NOT_POS_DEFINITE   -6

#define UQ_RC_MACRO(iRC,rank,where,what,retValue) \
  if (iRC) {                                      \
    std::cerr << "UQ RC ERROR"                    \
              << ", rank "  << rank               \
              << ", in "    << where              \
              << ": "       << what               \
              << ", iRC = " << iRC                \
              << std::endl;                       \
    return retValue;                              \
  }

#define UQ_TEST_MACRO(test,rank,where,what,retValue) \
  if (test) {                                        \
    std::cerr << "UQ TEST ERROR"                     \
              << ", rank "  << rank                  \
              << ", in " << where                    \
              << ": "    << what                     \
              << std::endl;                          \
    return retValue;                                 \
  }

#define UQ_FATAL_RC_MACRO(iRC,rank,where,what) \
  if (iRC) {                                   \
    std::cerr << "UQ RC FATAL ERROR"           \
              << ", rank "  << rank            \
              << ", in "    << where           \
              << ": "       << what            \
              << ", iRC = " << iRC             \
              << ". Exiting..."                \
              << std::endl;                    \
    exit(1);                                   \
  }

#define UQ_FATAL_TEST_MACRO(test,rank,where,what) \
  if (test) {                                     \
    std::cerr << "UQ TEST FATAL ERROR"            \
              << ", rank "  << rank               \
              << ", in "    << where              \
              << ": "       << what               \
              << ". Exiting..."                   \
              << std::endl;                       \
    exit(1);                                      \
  }

#endif // __UQ_DEFINES_H__
