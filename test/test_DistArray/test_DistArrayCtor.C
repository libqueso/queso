#include <stdio.h>
#include <string.h>
#include <queso/DistArray.h>

int main(int argc, char **argv) {
  QUESO::DistArray<std::string> d;

  /*
   * This code should never get here. If it does, the bash script that wraps
   * around it negates the return value, making this a failure
   */
  return 0;
}
