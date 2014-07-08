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
  std::size_t nlocalverts;//number of vertices stored on the local rank
  std::size_t nlocaledges;//number of edges stored on the local rank
  std::int64_t nglobalverts;//total number of vertices in the graph
  std::size_t *rowstarts;
  std::int64_t *column;
  //zero-based compressed sparse row data structure for the local part of graph
} csr_graph;

//void setup_globals(void); /* In utils.cpp */
void free_csr_graph(csr_graph* const g); /* In utils.cpp */
#endif
