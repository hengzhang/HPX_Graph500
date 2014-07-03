#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <mpi.h>
#include <assert.h>
#include "common.hpp"

//std::size_t rank, size;
void setup_globals(){
  //TODO: make multi localities - set rank and size of each locality

}

void free_csr_graph(csr_graph* const g){
  if (g->rowstarts != NULL) {free(g->rowstarts); g->rowstarts = NULL;}
  if (g->column != NULL) {free(g->column); g->column = NULL;}
}

