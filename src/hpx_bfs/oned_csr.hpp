#ifndef ONED_CSR_HPP
#define ONED_CSR_HPP

#include "common.hpp"
#include <limits.h>

#define ULONG_BITS (sizeof(unsigned long) * CHAR_BIT) 

typedef struct oned_csc_graph { 
  std::size_t nlocalverts; 
  int64_t max_nlocalverts; 
  std::size_t nlocaledges; 
  int lg_nglobalverts;
  int64_t nglobalverts; 
  std::size_t *rowstarts;
  int64_t *column;
  int lg_local_queue_size;
  const tuple_graph* tg;
} oned_csc_graph;
#define SWIZZLE_VERTEX(c) (((VERTEX_OWNER(c) <<lg_local_queue_size) * ULONG_BITS) | VERTEX_LOCAL(c))
#define UNSWIZZLE_VERTEX(c) (MUL_SIZE((c) & ((INT64_C(1)  << lg_local_queue_size) * ULONG_BITS - 1)) | ((((c) / ULONG_BITS) >> lg_local_queue_size)))

void convert_graph_to_oned_csc(const tuple_graph* const tg, oned_csc_graph* const g);

 void free_oned_csc_graph(oned_csc_graph* const g);

#endif
