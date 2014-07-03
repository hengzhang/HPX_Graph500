#ifndef COMMON_HPP
#define COMMON_HPP
#include <stdint.h>
#include <stddef.h>
#include <limits.h>

#include "../generator/graph_generator.hpp"

extern std::size_t rank , size;


#define MOD_SIZE(v) ((v) % size)
#define DIV_SIZE(v) ((v) / size)
#define MUL_SIZE(x) ((x) * size)
#define VERTEX_OWNER(v) ((int)(MOD_SIZE(v)))
#define VERTEX_LOCAL(v) ((size_t)(DIV_SIZE(v)))
#define VERTEX_TO_GLOBAL(i) ((int64_t)((i)*size + rank))

typedef struct csr_graph {
  size_t nlocalverts;
  size_t nlocaledges;
  int64_t nglobalverts;
  size_t *rowstarts;
  int64_t *column;
} csr_graph;

//void setup_globals(void); /* In utils.cpp */
void free_csr_graph(csr_graph* const g); /* In utils.cpp */
//void* xMPI_Alloc_mem(size_t nbytes); /* In utils.cpp */
//void* xmalloc(size_t nbytes); /* In utils.cpp */
//void* xcalloc(size_t n, size_t unit); /* In utils.cpp */
//void* xrealloc(void* p, size_t nbytes); /* In utils.cpp */

//void convert_graph_to_csr(const int64_t nedges, const int64_t* const edges, csr_graph* const g); /* In convert_to_csr.cpp */
//void find_bfs_roots(int *num_bfs_roots, const csr_graph* const g, const uint64_t seed1, const uint64_t seed2, int64_t* const bfs_roots); /* In find_roots.cpp */
//int validate_bfs_result(const csr_graph* const g, const int64_t root, const int64_t* const pred, const int64_t nvisited); /* In validate.cpp */

//void run_hpx_bfs(const csr_graph* const g, int64_t root, int64_t* pred, int64_t* nvisited); /* Provided by user */

#endif
