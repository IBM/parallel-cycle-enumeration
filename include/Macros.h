#ifndef _MACROS_H_
#define _MACROS_H_

#define MPI_IMPL

/** Define whether to use the improvements to the (temporal) Read-Tarjan algorithm **/
#define BLK_FORWARD
#define PATH_FORWARD
#define NEIGHBOR_ORDER
#define SUCCESSFUL_DFS_BLK

#ifndef BLK_FORWARD
#undef NEIGHBOR_ORDER
#undef SUCCESSFUL_DFS_BLK
#endif

#endif //_MACROS_H_
