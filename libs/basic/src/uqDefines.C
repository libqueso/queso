#include <uqDefines.h>

int uqMyRank() {
  int result = 0;
#ifdef __APPL_USES_TRILINOS__
  int iRC;
  iRC = MPI_Comm_rank(MPI_COMM_WORLD,&result);
#endif
  return result;
}
