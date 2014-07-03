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

void setup_globals(void); /* In utils.cpp */
void free_csr_graph(csr_graph* const g); /* In utils.cpp */
//void* xmalloc(size_t nbytes); /* In utils.cpp */
//void* xcalloc(size_t n, size_t unit); /* In utils.cpp */
//void* xrealloc(void* p, size_t nbytes); /* In utils.cpp */

#endif
