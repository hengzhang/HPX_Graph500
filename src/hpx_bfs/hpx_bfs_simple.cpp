#include "common.hpp"
#include "oned_csr.hpp"
#include <hpx/hpx.hpp>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>


/*static oned_csr_graph g;
static boost::uint64_t* g_oldq;
static boost::uint64_t* g_newq;
static unsigned long* g_visited;
static const int coalescing_size = 256;
static boost::uint64_t* g_outgoing;
static std::size_t* g_outgoing_counts 
static MPI_Request* g_outgoing_reqs;
static int* g_outgoing_reqs_active; 
static boost::uint64_t* g_recvbuf;
*/

//make this as a stub _ run bfs_custom
namespace graph500{ namespace bfs{

void find_bfs_roots(int* num_bfs_roots, const csr_graph* const g, const boost::uint64_t seed1, const uboost::uint64_t seed2, boost::uint64_t* const bfs_roots) {
  /* This implementation is slow, but there aren't enough roots being
   * generated for that to be a big issue. */
}


//!!!!just start run_bfs function
/*
* This version is the traditional level-synchronized BFS using two queues.  A
* bitmap is used to indicate which vertices have been  visited.  Messages are
* sent and processed asynchronously throughout the code  to hopefully overlap
* communication with computation.
*/
void run_hpx_bfs(const csr_graph *g, boost::uint64_t root, boost::uint64_t *pred, boost::uint64_t *nvisited){
  //implement in partition_component

}

}}
